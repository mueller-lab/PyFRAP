#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Geometry class for PyFRAP toolbox, including following classes:

#(1) geometry
#(2) zebrafishDomeStage
#(3) cylinder
#(4) xenopusBall

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_gmsh_module
from pyfrp.modules import pyfrp_gmsh_IO_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules.pyfrp_term_module import *

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================

class geometry(object):

	def __init__(self,embryo,typ,fnGeo,center):
		
		self.typ=typ
		self.embryo=embryo
		self.fnGeo=fnGeo
		self.center=center
		
		self.geoFileParameters={}
		
		
	def getEmbryo(self):
		return self.embryo
	
	def getTyp(self):
		return self.typ
	
	def setFnGeo(self,fn):
		self.fnGeo=fn
		return self.fnGeo
	
	def getFnGeo(self):
		return self.fnGeo
	
	def setCenter(self,c):
		self.center=c
		return self.center
	
	def getCenter(self):
		return self.center
	
	def center2Mid(self):
		if mod(self.embryo.dataResPx,2)==0:
			return self.setCenter([self.embryo.dataResPx/2+0.5,self.embryo.dataResPx/2+0.5])
		else:
			return self.setCenter([self.embryo.dataResPx/2,self.embryo.dataResPx/2])
	
	def centerInImg(self):
		oldCenter=self.getCenter()
		self.centerMid()
		a=raw_input("Change center of geometry from " + oldCenter + ' to ' + self.getCenter() + ' ? [Y/N]')
		if a=='Y':
			self.updateGeoFile()
			if None!=self.embryo.simulation:
				self.embryo.simulation.mesh.genMesh()
		
		else:
			self.setCenter(oldCenter)
		
		return 	self.getCenter()
	
	def setAllROI(self,name='All',makeNew=False,updateIdxs=False):
		
		if not hasattr(self,'optimalAllROI'):
			printWarning('Cannot set AllROI. Geometry of type ' + self.typ + ' has no method for this!')
			return None
		
		if makeNew:
			rnew=self.optimalAllROI(name)
		else:	
			r=self.embryo.getROIByName(name)
			if r==None:
				printWarning('Cannot set AllROI. ROI with name ' + name + ' does not exist!')
				return None
			
			
			incl=r.findIncluded()
				
			name=r.getName()
			Id=r.getId()
			color=r.getColor()
			asMaster=r.isMaster()
			
			rnew=self.optimalAllROI(name=name,Id=Id,color=color,asMaster=asMaster)
			
			for rincl in incl:
				ind=rincl.ROIsIncluded.index(r)
				rincl.ROIsIncluded.pop(ind)
				rincl.ROIsIncluded.insert(ind,rnew)
				
				if updateIdxs:
					rincl.updateIdxs()
			
			ind=self.embryo.ROIs.index(r)
			self.embryo.removeROI(ind)
			self.embryo.ROIs.insert(ind,self.embryo.ROIs.pop(-1))
				
		return rnew
	
	def readGeoFile(self):
		
		try:
			domain,self.geoFileParameters=pyfrp_gmsh_IO_module.readGeoFile(self.fnGeo)
		except:
			printWarning("Cannot read geo file " + self.fnGeo)
			return None
		
		return domain
	
	def plotGeometry(self,ax=None,color='k',ann=False):
		
		domain=self.readGeoFile()
		
		if domain==None:
			return ax
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Geometry" +self.typ])
			ax=axes[0]
		
		domain.draw(ax=ax,color=color,ann=ann)
		
		return ax
		
class zebrafishDomeStage(geometry):
	
	def __init__(self,embryo,center,imagingRadius,radiusScale=1.1):
		
		super(zebrafishDomeStage, self).__init__(embryo,"zebrafishDomeStage","../meshfiles/dome.geo",center)
		
		#How much bigger is the inner radius than the outerRadius
		self.radiusScale=radiusScale
		self.imagingRadius=imagingRadius
	
		#Compute geometry properties
		self.restoreDefault()
	
	def setOuterRadius(self,r):
		self.outerRadius=r
		self.computeDome()
		return self.outerRadius
	
	def setRadiusScale(self,s):
		self.radiusScale=s
		self.computeDome()
		return self.radiusScale
	
	def setImagingRadius(self,r):
		self.imagingRadius=r
		self.computeDome()
		return self.imagingRadius
	
	def setImagingHeight(self,h):
		self.imagingHeight=h
		if self.embryo.sliceHeightPx!=self.imagingHeight:
			printWarning("imagingHeight of geometry is not identical with imaging depth of embryo. This can lead to problems.")
		
		self.computeDome()
		return self.imagingHeight
	
	def getImagingRadius(self):
		return self.imagingRadius
	
	def getImagingHeight(self):
		return self.imagingHeight
	
	def getOuterRadius(self):
		return self.outerRadius
	
	def getRadiusScale(self):
		return self.radiusScale
	
	def computeDome(self):
		self.outerRadius=(self.imagingRadius**2+self.imagingHeight**2)/(2*(-self.imagingHeight))
		self.innerRadius=self.outerRadius*self.radiusScale
		self.centerDist=np.sqrt(self.innerRadius**2-self.outerRadius**2) 
		return self.innerRadius,self.centerDist
	
	def restoreDefault(self):
		self.imagingHeight=self.embryo.sliceHeightPx
		self.computeDome()
	
	def updateGeoFile(self,debug=False):
		pyfrp_gmsh_module.updateDomeGeo(self.fnGeo,self.imagingRadius,self.imagingHeight,self.center,run=False,debug=debug)
		
	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getOuterRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
		
class zebrafishDomeStageQuad(zebrafishDomeStage):
	
	def __init__(self,embryo,center,imagingRadius,radiusScale=1.1):
		
		zebrafishDomeStage.__init__(self,embryo,center,imagingRadius,radiusScale=1.1)
		geometry.__init__(self,embryo,"zebrafishDomeStageQuad","../meshfiles/quad_dome.geo",center)
	
	
class cylinder(geometry):
	
	def __init__(self,embryo,center,radius,height):
		
		super(cylinder, self).__init__(embryo,"cylinder","../meshfiles/cylinder.geo",center)

		self.radius=radius
		self.height=height
		
	def setHeight(self,h):
		self.height=h
		return self.height
	
	def setRadius(self,r):
		self.radius=r
		return self.radius
	
	def getHeight(self):
		return self.height
	
	def getRadius(self):
		return self.radius
	
	def updateGeoFile(self,debug=False):
		pyfrp_gmsh_module.updateCylinderGeo(self.fnGeo,self.radius,self.height,self.center,run=False,debug=debug)

	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
	
class cylinderQuad(cylinder):
	
	def __init__(self,embryo,center,radius,height):
		
		cylinder.__init__(self,embryo,center,radius,height)
		geometry.__init__(self,embryo,"cylinderQuad","../meshfiles/quad_cylinder.geo",center)
		
	
class xenopusBall(geometry):
	
	def __init__(self,embryo,center,imagingRadius):		
		
		super(xenopusBall, self).__init__(embryo,"xenopusBall","../meshfiles/ball.geo",center)
		
		self.imagingRadius=imagingRadius
		
		self.restoreDefault()
	
	def setImagingRadius(self,r):
		self.imagingRadius=r
		self.computeDome()
		return self.imagingRadius
	
	def setImagingHeight(self,h):
		self.imagingHeight=h
		self.computeDome()
		return self.imagingHeight
	
	def getImagingRadius(self,r):
		return self.imagingRadius
	
	def getImagingHeight(self,h):
		return self.imagingHeight
	
	def computeBall(self):
		self.radius=(self.imagingRadius**2+self.imagingHeight**2)/(2*(-self.imagingHeight))
	
	def restoreDefault(self):
		self.imagingHeight=self.embryo.sliceHeightPx
		self.computeBall()
		
	def updateGeoFile(self,debug=False):
		pyfrp_gmsh_module.updateBallGeo(self.fnGeo,self.radius,self.center,run=False,debug=debug)
	
	def getRadius(self):
		return self.radius
	
	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
	
class xenopusBallQuad(xenopusBall):
	
	def __init__(self,embryo,center,imagingRadius):
		
		xenopusBall.__init__(self,embryo,center,imagingRadius)	
		geometry.__init__(self,embryo,"xenopusBallQuad","meshfiles/quad_ball.geo",center)	

class cone(geometry):
	
	def __init__(self,embryo,center,upperRadius,lowerRadius,height):
		
		super(cone, self).__init__(embryo,"cone","../meshfiles/cone.geo",center)

		self.upperRadius=upperRadius
		self.lowerRadius=lowerRadius
		
		self.height=height
		
	def setHeight(self,h):
		self.height=h
		return self.height
	
	def setLowerRadius(self,r):
		self.lowerRadius=r
		return self.lowerRadius
	
	def getLowerRadius(self):
		return self.lowerRadius
	
	def setUpperRadius(self,r):
		self.upperRadius=r
		return self.upperRadius
	
	def getUpperRadius(self):
		return self.upperRadius
	
	def getHeight(self):
		return self.height
	
	def computeSliceHeightFromRadius(self,radius):
		sliceHeight=self.height/(self.upperRadius-self.lowerRadius)*(radius-self.lowerRadius) 
		return sliceHeight
	
	def computeRadiusFromSliceHeight(self,height):
		radius=(self.upperRadius-self.lowerRadius)/self.height*height+self.upperRadius
		return radius
	
	def updateGeoFile(self,debug=False):
		pyfrp_gmsh_module.updateConeGeo(self.fnGeo,self.upperRadius,self.lowerRadius,abs(self.height),self.center,run=False,debug=debug)

	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getUpperRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
		
class custom(geometry):

	def __init__(self,embryo,center,fnGeo):
		
		geometry.__init__(self,embryo,"custom",fnGeo,center)
	
	