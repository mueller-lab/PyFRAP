#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFRAP.

#PyFRAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#ROI module for PyFRAP toolbox, including following ROI objects:

#(1)  ROI
#(2)  radialROI
#(3)  sliceROI
#(4)  radialSliceROI
#(5)  squareROI
#(6)  squareSliceROI
#(7)  rectangleROI
#(8)  rectangleSliceROI
#(9)  polyROI
#(10)  polySliceROI
#(11)  customROI

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module 
from pyfrp.modules import pyfrp_idx_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_img_module
from pyfrp.modules import pyfrp_integration_module
from pyfrp.modules import pyfrp_fit_module

from pyfrp.modules.pyfrp_term_module import *

#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Time 
import time

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main ROI class

class ROI(object):

	def __init__(self,embryo,name,Id,zmin='-inf',zmax='inf',color='b'):
		
		#Name/Id
		self.name=name
		self.Id=Id
		self.embryo=embryo
		self.color=color
		
		#zExtend
		self.zmin=pyfrp_misc_module.translateNPFloat(zmin)
		self.zmax=pyfrp_misc_module.translateNPFloat(zmax)
		
		#Idxs from data analysis/simulation
		self.imgIdxX=[]
		self.imgIdxY=[]
		self.extImgIdxX=[]
		self.extImgIdxY=[]
		self.meshIdx=[]
		
		#Mask
		self.imgMask=None
		self.extMask=None
		
		#Number of extended pixels
		self.numExt=None
		
		#Result Dataseries
		self.dataVec=[]
		self.simVec=[]
		self.dataVecPinned=[]
		self.simVecPinned=[]
		
		#Rim concentration
		self.useForRim=False
	
	def setName(self,n):
		self.name=n
		return self.name
	
	def setZExtend(self,zmin,zmax):
		
		self.zmin=zmin
		self.zmax=zmax
		
		return [self.zmin,self.zmax]
	
	def getZExtend(self):
		return [self.zmin,self.zmax]
	
	def getId(self):
		return self.Id
	
	def getName(self):
		return self.name
	
	def getImgIdx(self):
		return self.imgIdxX,self.imgIdxY
	
	def getExtImgIdx(self):
		return self.extImgIdxX,self.extImgIdxY
	
	def getMeshIdx(self):
		return self.meshIdx
	
	def getMeshDensity(self):
		volume=self.getVolume()
		return len(self.meshIdx)/volume
		
	def getVolume(self):
		area=self.getArea()
		height=self.getROIHeight()
		
		return area*height
		
	def getROIHeight(self):
	
		if np.isfinite(self.zmax):
			zMax=self.zmax
		else:
			dump,zMax=self.getMeshIdxZExtend()
		
		if np.isfinite(self.zmin):
			zMin=self.zmin
		else:
			zMin,dump=self.getMeshIdxZExtend()
			
		return abs(zMax-zMin)
	
	def getArea(self):
		if self.imgMask==None:
			self.computeImgMask()
		if self.extMask==None:
			self.computeExtMask()
			if self.extMask==None:
				return self.imgMask.sum()
		
		return self.imgMask.sum()+self.extMask.sum()
		
	def getMeshIdxZExtend(self):
		
		"""Returns extend of ROI's ``meshIdx`` in z-coordinate.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum z-coordinate.
				* (float): Maximum z-coordinate.
				
		"""
		
		mesh=self.embryo.simulation.mesh.mesh
		z=np.asarray(mesh.z)[self.meshIdx]	
		
		return min(z) , max(z)
	
	def getMeshIdxYExtend(self):
		
		"""Returns extend of ROI's ``meshIdx`` in y-coordinate.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum y-coordinate.
				* (float): Maximum y-coordinate.
				
		"""
		
		mesh=self.embryo.simulation.mesh.mesh
		y=np.asarray(mesh.y)[self.meshIdx]	
		
		return min(y) , max(y)
	
	def getMeshIdxXExtend(self):
		
		"""Returns extend of ROI's ``meshIdx`` in x-coordinate.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum x-coordinate.
				* (float): Maximum x-coordinate.
				
		"""
		
		mesh=self.embryo.simulation.mesh.mesh
		x=np.asarray(mesh.x)[self.meshIdx]	
		
		return min(x) , max(x)
	
	def getMeshIdxExtend(self):
		
		"""Returns extend of ROI's ``meshIdx``.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum x-coordinate.
				* (float): Maximum x-coordinate.
				* (float): Minimum y-coordinate.
				* (float): Maximum y-coordinate.
				* (float): Minimum z-coordinate.
				* (float): Maximum z-coordinate.
		"""
		
		xmin,xmax=self.getMeshIdxXExtend()
		ymin,ymax=self.getMeshIdxYExtend()
		zmin,zmax=self.getMeshIdxZExtend()
		
		return xmin,xmax,ymin,ymax,zmin,zmax
		
	def getType(self):
		typ=str(type(self))
		before,typ,after=typ.split("'")
		typ=typ.replace('pyfrp_ROI.','')
		typ=typ.replace('ROI','')
		
		typ=typ.replace('pyfrp.subclasses.','')
		
		return typ
	
	def setColor(self,color):
		self.color=color
		return self.color
	
	def setUseForRim(self,b):
		if self.numExt>0 and b==True:
			printWarning('Be careful, region '+self.name+' is set for rim calculation but has indices outside of image.')
		self.useForRim=b
		return self.useForRim
	
	def getUseForRim(self):
		return self.useForRim
	
	def getColor(self):
		return self.color
	
	def emptyIdxs(self):
		
		self.imgIdxX=[]
		self.imgIdxY=[]
		self.meshIdx=[]
		
		return self.getAllIdxs()
	
	def copyIdxs(self,r):
		self.imgIdxX=r.imgIdxX
		self.imgIdxY=r.imgIdxY
		
		self.extImgIdxX=r.extImgIdxX
		self.extImgIdxY=r.extImgIdxY
		
		self.meshIdx=r.meshIdx
		
		return self.getAllIdxs()
	
	def getAllIdxs(self):
		return self.imgIdxX,self.imgIdxY,self.meshIdx
	
	def getImgMask(self):
		return self.imgMask
	
	def getExtMask(self):
		return self.extMask
	
	def computeNumExt(self):
		self.numExt=len(self.extImgIdxX)
		return self.numExt
	
	def getNumExt(self):
		return self.numExt
	
	def setDataVec(self,vec):
		self.dataVec=vec 
		return self.dataVec

	def getDataVec(self):
		return self.dataVec
	
	def setSimVec(self,vec):
		self.simVec=vec 
		return self.simVec

	def getSimVec(self):
		return self.simVec
	
	def computeImgMask(self):
		vals=np.zeros((self.embryo.dataResPx,self.embryo.dataResPx))
		self.imgMask=pyfrp_idx_module.ind2mask(vals,self.imgIdxX,self.imgIdxY,1)
		return self.imgMask
		
	def computeExtMask(self):
		
		if len(self.extImgIdxX)==0:
			return None,None,None
		
		minX=min(self.extImgIdxX)
		maxX=max(self.extImgIdxX)
		
		minY=min(self.extImgIdxY)
		maxY=max(self.extImgIdxY)
		
		X=np.arange(minX,maxX+1)
		Y=np.arange(minY,maxY+1)
		
		mX,mY=np.meshgrid(X,Y)
		
		vals=np.zeros((len(X),len(Y)))
		
		idXtemp=np.asarray(self.extImgIdxX)+abs(minX)
		idYtemp=np.asarray(self.extImgIdxY)+abs(minY)
		
		idXtemp=idXtemp.astype('int')
		idYtemp=idYtemp.astype('int')
		
		self.extMask=pyfrp_idx_module.ind2mask(vals,idXtemp,idYtemp,1)
		
		return mX,mY,self.extMask
		
	def showImgIdx(self,ax=None):
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["imgIdx"],sup=self.name+" imgIdx")
			ax=axes[0]
		
		self.computeImgMask()
	
		ax.imshow(self.imgMask)
		plt.draw()
		
		
		return ax
	
	def showExtImgIdx(self,ax=None):
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["extImgIdx"],sup=self.name+" imgIdx")
			ax=axes[0]
		
		mX,mY,self.extMask=self.computeExtMask()
		if mX!=None:
			ax.contourf(mX,mY,self.extMask)
			plt.draw()
			
		return ax
	
	def showMeshIdx(self,ax=None):
		
		mesh=self.embryo.simulation.mesh.mesh
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["MeshIdx"],sup=self.name+" MeshIdx",proj=['3d'])
			ax=axes[0]
		
		#Somehow need to convert to np array since slicing does not work for fipy variables
		x=np.asarray(mesh.x)[self.meshIdx]
		y=np.asarray(mesh.y)[self.meshIdx]
		z=np.asarray(mesh.z)[self.meshIdx]
		
		ax.scatter(x,y,z,c=self.color)
		plt.draw()
		
		return ax
	
	def showMeshIdx2D(self,ax=None):
		
		mesh=self.embryo.simulation.mesh.mesh
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["MeshIdx"],sup=self.name+" MeshIdx")
			ax=axes[0]
		
		#Somehow need to convert to np array since slicing does not work for fipy variables
		x=np.asarray(mesh.x)[self.meshIdx]
		y=np.asarray(mesh.y)[self.meshIdx]
	
		ax.scatter(x,y,c=self.color)
		plt.draw()
		
		return ax
	
	
	def computeIdxs(self,matchMesh=False,debug=False):
		
		"""Computes image and mesh indices of ROI. 
	
		Will do this by:
		
			* Compute image indices.
			* Match image indices with master ROI.
			* Compute external indices.
			* Compute mesh indices.
			* Match mesh indices with the ones of master ROI.
			
		.. note:: If no master ROI is defined, will not do anything.
		
		.. note:: If master ROI has not been indexed yet, will first index it, then continue. 
		
		.. note:: Will skip mesh index computation if there is no mesh generated yet.
		
		Keyword Args:
			matchMesh (bool): Match mesh indices with master ROI.
			debug (bool): Print out debugging messages.
			
		Return:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x direction.
				* imgIdxY (list): Image indices in y direction.
				* meshIdx (list): Mesh indices.
		
		"""
		
		if self.embryo.getMasterROIIdx()==None:
			printWarning("No Master ROI has been defined yet. Will not continue compute ROI indices.")
			return self.getAllIdxs()
		else:
			masterROI=self.embryo.getMasterROI()
		
		if self!=masterROI:
			if len(masterROI.imgIdxX)==0:
				printWarning("Idxs of Master ROI have not been computed. Will compute them first.")
				masterROI.computeIdxs(debug=debug)
		
		startIdx=time.clock()
		
		if  type(self) is not customROI:
			
			self.computeImgIdx(debug=debug)
			self.matchImgIdx(masterROI)
			self.computeExtIdx(debug=debug)
			
			if self.embryo.simulation!=None:
				if self.embryo.simulation.mesh.mesh==None:
					printWarning("Mesh has not been generated, will not compute meshIdxs")
				else:
				
					self.computeMeshIdx(self.embryo.simulation.mesh.mesh)
				
					if matchMesh:
						if self!=masterROI:
							self.matchMeshIdx(masterROI)
			else:
				printWarning("Simulation object does not exist yet, hence won't index for mesh.")
		else:
			
			self.updateIdxs()
			self.matchImgIdx(masterROI)
			if self.embryo.simulation!=None:
				if self.embryo.simulation.mesh.mesh!=None:
					self.matchMeshIdx(masterROI)
		
		if debug:
			print 'Compute Idxs: ', startIdx-time.clock()
		
		return self.getAllIdxs()
	
	def computeExtIdx(self,debug=False):
		
		"""Computes indices of external pixels.
		
		Does this by comparing extended pixels of ``self`` with the one of the master ROI.
		
		Keyword Args:
			debug (bool): Print out debugging messages.
		
		Return:
			tuple: Tuple containing:
			
				* extImgIdxX (list): External image indices in x direction.
				* extImgIdxY (list): External image indices in y direction.
		
		"""
		
		m=self.embryo.getMasterROI()
		rois=[self,m]
		[self.extImgIdxX,self.extImgIdxY]=pyfrp_idx_module.getCommonExtendedPixels(rois,self.embryo.dataResPx,debug=debug)
		return self.extImgIdxX,self.extImgIdxY
	
	def matchImgIdx(self,r):
		
		"""Matches image indices of ``self`` with the ones of ROI ``r``.
		
		Does this by generating masks of both ROIs and multiplicating them.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to match with.
		
		Return:
			tuple: Tuple containing:
			
				* imgIdxX (list): Matched image indices in x direction.
				* imgIdxY (list): Matched image indices in y direction.
		
		"""
		
		self.computeImgMask()
		self.imgMask=self.imgMask*r.computeImgMask()
		self.imgIdxX,self.imgIdxY=pyfrp_idx_module.mask2ind(self.imgMask,self.embryo.dataResPx)	
		return self.imgIdxX,self.imgIdxY
	
	def matchMeshIdx(self,r,matchZ=False):
		
		x=np.asarray(self.embryo.simulation.mesh.mesh.x)[self.meshIdx]
		y=np.asarray(self.embryo.simulation.mesh.mesh.y)[self.meshIdx]
		
		ins=r.checkXYInside(x,y)
		
		self.meshIdx=np.asarray(self.meshIdx)
		self.meshIdx=self.meshIdx[np.where(ins)[0]]
		
		return self.meshIdx
	
	def checkZInside(self,z):
		if self.zmin<= z and z<=self.zmax:
			return True
		else:
			return False
	
	def showIdxs(self,axes=None):
		wereNone=False
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["MeshIdx","imgIdx","extImgIdx"],sup=self.name+" Idx",proj=['3d',None,None])
			wereNone=True
			
		self.showMeshIdx(axes[0])
		self.showImgIdx(axes[1])
		self.showExtImgIdx(axes[2])
		
		if wereNone:
			axes[0].set_title(self.name+" Idx")
			
		return axes
	
	def checkSymmetry(self,debug=False):
		img=np.zeros((self.embryo.dataResPx,self.embryo.dataResPx))
		img[self.imgIdxX,self.imgIdxY]=1
		return pyfrp_img_module.symmetryTest(img,debug=debug)
		
	def idxs2Quad(self,debug=False):
		
		if not self.checkSymmetry():
			printWarning('Cannot reduce region '+self.name+' to quadrant. Indices are not symmetric.')
			return self.getAllIdxs()
		
		self.imgIdxX,self.imgIdxY=imgIdx2Quad(self.imgIdxX,self.imgIdxY,self.embryo.dataResPx,debug=debug)
		
		if 'Quad' not in self.embryo.geometry.typ:
			printWarning('Will not update mesh indices, geometry is not set to quad.')
		else:	
			self.computeMeshIdx()
		
		return self.getAllIdxs()
	
	def idxs2Full(self):
		return self.computeIdxs()
	
	def resetDataVec(self):
		self.setDataVec([])
		return self
	
	def resetSimVec(self):
		self.setSimVec([])
		return self
	
	def plotData(self,ax=None,color=None,linewidth=1,legend=True,linestyle='-'):
		
		if color==None:
			color=self.color
		
		ax = pyfrp_plot_module.plotTS(self.embryo.tvecData,self.dataVec,ax=ax,linewidth=linewidth,color=color,label=self.name + ' data',
		title="Data",sup=self.name+" data",linestyle=linestyle,legend=legend)
		
		return ax
	
	def plotDataPinned(self,ax=None,color=None,linewidth=1,legend=True,linestyle='-'):
		
		if color==None:
			color=self.color
		
		ax = pyfrp_plot_module.plotTS(self.embryo.tvecData,self.dataVecPinned,ax=ax,linewidth=linewidth,color=color,label=self.name + ' data',
		title="Data Pinned",sup=self.name+" data",linestyle=linestyle,legend=legend)
		
		return ax
	
	def plotSim(self,ax=None,color=None,linewidth=1,legend=True,linestyle='--'):
		
		if color==None:
			color=self.color
		
		ax = pyfrp_plot_module.plotTS(self.embryo.simulation.tvecSim,self.simVec,ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' simulated',title="Simulation",sup=self.name+" simulation",linestyle=linestyle,legend=legend)
			
		return ax
	
	def plotSimPinned(self,ax=None,color=None,linewidth=1,legend=True,linestyle='--'):
		
		if color==None:
			color=self.color
		
		ax = pyfrp_plot_module.plotTS(self.embryo.simulation.tvecSim,self.simVecPinned,ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' ' + ' simulated',title="Simulation Pinned",sup=self.name+" simulation",linestyle=linestyle,legend=legend)
			
		return ax

	def findIncluded(self):
		incl=[]
		for r in self.embryo.ROIs:
			if type(r) is customROI:
				if r.roiIncluded(self):
					incl.append(r)
		return incl			
					
	def isMaster(self):
		return self==self.embryo.getMasterROI()
	
	def plotSolutionVariable(self,phi,ax=None,vmin=None,vmax=None,nlevels=25,colorbar=True):
		
		if hasattr(phi,'value'):
			val=phi.value
		else:
			val=phi
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Solution Variable"],sup=self.name+" phi")
			ax=axes[0]
		
		x,y,z=self.embryo.simulation.mesh.mesh.getCellCenters()
		
		if vmin==None:
			vmin=min(val[self.meshIdx])
		if vmax==None:
			vmax=max(val[self.meshIdx])
		
		levels=np.linspace(vmin,1.01*vmax,nlevels)
		
		solPlot=ax.tricontourf(x[self.meshIdx],y[self.meshIdx],val[self.meshIdx],vmin=vmin,vmax=vmax,levels=levels)
		ax.autoscale(enable=True, axis='both', tight=True)
		
		if colorbar:
			cb=plt.colorbar(solPlot,orientation='horizontal',pad=0.05)
		
		ax.get_figure().canvas.draw()
		
		return ax
	
	def getSimConc(self,phi,append=True):
		
		cvs=self.embryo.simulation.mesh.mesh.getCellVolumes()
		
		c=pyfrp_integration_module.getAvgConc(phi,cvs,self.meshIdx)
		
		if append:
			self.simVec.append(c)
			
		return c	
	
	def pinAllTS(self,bkgdVal=None,normVal=None,bkgdValSim=None,normValSim=None,debug=False):
		
		bkgdValSim=pyfrp_misc_module.assignIfVal(bkgdValSim,bkgdVal,None)
		normValSim=pyfrp_misc_module.assignIfVal(normValSim,normVal,None)
		
		self.dataVecPinned = self.pinDataTS(bkgdVal=bkgdVal,normVal=normVal,debug=debug)
		self.simVecPinned = self.pinSimTS(bkgdVal=bkgdValSim,normVal=normValSim,debug=debug)
		
		return self.dataVecPinned,self.simVecPinned
	
	def pinDataTS(self,bkgdVal=None,normVal=None,debug=False):
		
		if bkgdVal==None or normVal==None:
			bkgdValTemp,normValTemp = self.embryo.computePinVals(debug=debug)
			bkgdVal = pyfrp_misc_module.assignIfVal(bkgdVal,bkgdValTemp,None)
			normVal = pyfrp_misc_module.assignIfVal(normVal,normValTemp,None)
			
		self.dataVecPinned=pyfrp_fit_module.pinConc(self.dataVec,bkgdVal,normVal,axes=None,debug=debug,tvec=self.embryo.tvecData,color=self.color)
		
		return self.dataVecPinned
	
	def pinSimTS(self,bkgdVal=None,normVal=None,debug=False):
		
		if bkgdVal==None or normVal==None:
			bkgdValTemp,normValTemp = self.embryo.computePinVals(debug=debug)
			bkgdVal = pyfrp_misc_module.assignIfVal(bkgdVal,bkgdValTemp,None)
			normVal = pyfrp_misc_module.assignIfVal(normVal,normValTemp,None)
		
		self.simVecPinned=pyfrp_fit_module.pinConc(self.simVec,bkgdVal,normVal,axes=None,debug=debug,tvec=self.embryo.simulation.tvecSim,color=self.color)
		
		return self.simVecPinned
	
	def getFittedVec(self,fit):
		try:
			return fit.fittedVecs[fit.ROIsFitted.index(self)]
		except IndexError:
			fit.fittedVecs.insert(fit.ROIsFitted.index(self),[])
			#return fit.fittedVecs[fit.ROIsFitted.index(self)]	
			#This solves a problem only temporarily.
			return []
		
	def getdataVecFitted(self,fit):
		return fit.dataVecsFitted[fit.ROIsFitted.index(self)]
		
	def plotFit(self,fit,ax=None,color=None,linewidth=1,legend=True,title=None):
		
		if color==None:
			color=self.color
		
		if title==None:
			title="Fit"+fit.name
			
		ax = pyfrp_plot_module.plotTS(fit.tvecFit,self.getFittedVec(fit),ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' ' + fit.name,title=title,sup=self.name+" fitted",linestyle='-.',legend=legend)
		
		ax = pyfrp_plot_module.plotTS(fit.tvecFit,self.getdataVecFitted(fit),ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' ' + fit.name,title=title,sup=self.name+" fitted",linestyle='-',legend=legend)
		
		return ax
	
	def isAnalyzed(self):
		return len(self.embryo.tvecData)==len(self.dataVec)
	
	def isSimulated(self):
		if self.embryo.simulation!=None:
			return len(self.embryo.simulation.tvecSim)==len(self.simVec)
		return False
	
	def isFitted(self):
		fitted=False
		for fit in self.embryo.fits:
			if self in fit.ROIsFitted and len(self.getFittedVec(fit))>0:
				fitted=True
		return fitted
	
	def getInterpolationError(self):
		if self.isSimulated() and self.isAnalyzed():
			try:
				return self.dataVec[0]/self.simVec[0]
			except ZeroDivisionError:
				printWarning("Dividing by zero. Going to return infinity.")
				return np.inf
		else:
			printWarning("ROI is either not simulated or analyzed. Cannot return interpolation error.")
			return 0.
	
	def getEncapsulatingBox(self):
		xExtend,yExtend=self.computeXYExtend()
		zExtend=self.getZExtend()
		return xExtend,yExtend,zExtend
		
	def refineInMesh(self,factor=3.,addZ=15.,findIdxs=True,debug=False,run=True):
		
		xExtend,yExtend,zExtend=self.getEncapsulatingBox()
		zExtend=[zExtend[0]-addZ,zExtend[1]+addZ]
		
		print "current zExtend", zExtend
		
		
		if debug:
			print "Adding Box Field for ROI " + self.name
			print "Mesh Nodes in ROI before: ", len(self.meshIdx)
			
		fnOut=self.embryo.simulation.mesh.addBoxField(self.embryo.simulation.mesh.volSizePx/factor,xExtend,yExtend,zExtend,comment=self.name+" field",run=run)
	
		if findIdxs:
			self.computeMeshIdx(self.embryo.simulation.mesh.mesh)
		
		if debug and findIdxs:
			print "Mesh Nodes in ROI after: ", len(self.meshIdx)
			
		return fnOut
	
	def adaptRefineInMesh(self,nNodesReq,factor=3.,addZ=15.,debug=False):
		
		"""Refines mesh inside ROI adaptively until a given number of nodes inside ROI 
		is reached.
		
		Does this by:
			
			* Refining through :py:func:`refineInMesh`.
			* Computing mesh indices via :py:func:`computeMeshIdx`.
			* If number of nodes did not change, increase ``addZ``, else increase ``factor``.
			* Check if desired number of nodes is reached or not, if not, repeat.
		
		Args:
			nNodesReq (int): Desired number of nodes inside ROI.
		
		Keyword Args:
			factor (float): Refinement factor.
			addZ (float): Number of pixels added above and below ROI for box field.
			debug (bool): Print debugging messages.
			
		Returns:
			int: Final number of nodes in ROI.
			
		"""
		
		
		nNodes=len(self.meshIdx)
		nNodesAll=self.embryo.simulation.mesh.getNNodes()
		i=0
		while nNodes<nNodesReq:
			
			self.refineInMesh(factor=factor,addZ=addZ,findIdxs=True,debug=debug,run=True)
			
			nNodesNew=len(self.meshIdx)
			nNodesAllNew=self.embryo.simulation.mesh.getNNodes()
			
			if debug:
				print "Iteration ", i, ". "
				print "Total mesh nodes: ", nNodesAllNew
				print "Mesh Nodes in ROI before refinement: " , nNodes, " and after ", nNodesNew, "."
			
			if nNodesNew<nNodesReq:
				if nNodesAllNew==nNodesAll:
					if debug: 
						print "nNodesAll did not change, will increase addZ by 1. \n"
					addZ=addZ+1
					
				else:
					if debug:
						print "nNodes not large enough yet, will increase factor by 1. \n"
					factor=factor+1
			i=i+1
			
			nNodes=nNodesNew
			nNodesAll=nNodesAllNew
		return nNodes
				
	def printDetails(self):
		
		"""Prints out all attributes of ROI object."""
		
		print "ROI ", self.name, " details:"
		printAllObjAttr(self)
		print 
	
	def plotSimConcProfile(self,phi,ax=None,direction='x'):
		
		"""Plots concentration profile of solution variable in 
		single direction.
		
		Args:
			phi (fipy.CellVariable): Solution variable
			
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			direction (str): Direction to be plotted (x/y/z).
		
		Returns:
			matplotlib.axes: Matplotlib axes used for plotting.
		
		"""
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Concentration profile"],sup=self.name+" phi")
			ax=axes[0]
		
		if direction=='x':
			x=self.embryo.simulation.mesh.mesh.x
		elif direction=='y':
			x=self.embryo.simulation.mesh.mesh.y
		elif direction=='z':
			x=self.embryo.simulation.mesh.mesh.z
		else:
			printError('Direction '+ direction+ 'unknown. Will not plot.')
			return ax
		
		x=np.asarray(x)[self.meshIdx]
		v=np.asarray(phi.value)[self.meshIdx]
		
		xSorted,vSorted=pyfrp_misc_module.sortListsWithKey(v,x)
		
		pyfrp_plot_module.plotTS(xSorted,vSorted,color=self.color,label=self.name,legend=True)
		ax.set_xlabel(direction)
		ax.set_ylabel("Concentration")
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Radial ROI class

class radialROI(ROI):

	def __init__(self,embryo,name,Id,center,radius,color='b'):
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.radius=radius
		self.center=center
		
	def setRadius(self,r):
		self.radius=r
		return self.radius	
		
	def getRadius(self):
		return self.radius
	
	def setCenter(self,c):
		self.center=c
		return self.center
	
	def getCenter(self):
		return self.center	
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getCircleIdxImg(self.center,self.radius,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getCircleIdxMesh(self.center,self.radius,radius,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx	
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
			
		patch = ptc.Circle(self.center,self.radius,fill=False,linewidth=linewidth,color=color)
		ax.add_patch(patch)
		return ax
	
	def center2Mid(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return self.setCenter([self.embryo.dataResPx/2+0.5,self.embryo.dataResPx/2+0.5])
		else:
			return self.setCenter([self.embryo.dataResPx/2,self.embryo.dataResPx/2])
	
	def makeReducable(self,auto=False,debug=False):
		
		oldCenter=self.getCenter()
		self.center2Mid()
		
		if not auto:
			a=raw_input("Change center of region "+ self.name + " from " + str(oldCenter) + ' to ' + str(self.getCenter()) + ' ? [Y/N]')
		
			if a=='N':
				self.setCenter(oldCenter)
				return False
			elif a=='Y':
				pass
			
			if not self.checkSymmetry():
				printWarning('Cannot make region '+self.name+' reducable.')
				self.setCenter(oldCenter)
				return False
			return True
		return False
		
		
	def checkCentered(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return bool((self.center[0]==self.embryo.dataResPx/2+0.5) and (self.center[1]==self.embryo.dataResPx/2+0.5))
		else:
			return bool((self.center[0]==self.embryo.dataResPx/2.) and (self.center[1]==self.embryo.dataResPx/2.))
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsideCircle(x,y,self.center,self.radius)
	
	def computeXYExtend(self):
		self.xExtend=[self.center[0]-self.radius,self.center[0]+self.radius]
		self.yExtend=[self.center[1]-self.radius,self.center[1]+self.radius]
		return self.xExtend, self.yExtend
	

	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#slice ROI class
	
class sliceROI(ROI):
	
	def __init__(self,embryo,name,Id,height,width,sliceBottom,color='b'):
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.height=height
		self.width=width
		self.sliceBottom=sliceBottom
		
		self.computeZExtend()
		
	def computeZExtend(self):
		
		if self.sliceBottom:
			self.setZExtend(self.height,self.height+self.width)
		else:
			self.setZExtend(self.height-0.5*self.width,self.height+0.5*self.width)
		
		return self.zmin,self.zmax
		
	def setHeight(self,h):
		self.height=h
		self.computeZExtend()
		return self.height
	
	def getHeight(self):
		return self.height
	
	def setSliceBottom(self,s):
		self.sliceBottom=s
		self.computeZExtend()
		return self.sliceBottom
	
	def getSliceBottom(self):
		return self.sliceBottom
	
	def setWidth(self,w):
		self.width=w
		self.computeZExtend()
		return self.width
	
	def getWidth(self):
		return self.width
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getAllIdxImg(self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		x,y,z=mesh.cellCenters
		self.meshIdx=pyfrp_idx_module.getSliceIdxMesh(z,self.zmin,self.zmax)
		return self.meshIdx
	
	def checkXYInside(self,x,y):
		return True
	
	def computeXYExtend(self):
		self.xExtend=[0,self.embryo.dataResPx]
		self.yExtend=[0,self.embryo.dataResPx]
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Radial and slice ROI class
	
class radialSliceROI(sliceROI,radialROI):
	
	def __init__(self,embryo,name,Id,center,radius,height,width,sliceBottom,color='b'):
		
		radialROI.__init__(self,embryo,name,Id,center,radius,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
		
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getCircleIdxImg(self.center,self.radius,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
	
		self.meshIdx=pyfrp_idx_module.getCircleIdxMesh(self.center,self.radius,mesh,zmin=self.zmin,zmax=self.zmax)
		
		return self.meshIdx
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsideCircle(x,y,self.center,self.radius)
	
	def computeXYExtend(self):
		self.xExtend=[self.center[0]-self.radius,self.center[0]+self.radius]
		self.yExtend=[self.center[1]-self.radius,self.center[1]+self.radius]
		return self.xExtend, self.yExtend
	
	#def plotIn3D(self,domain=None,ax=None):
	
	###NOTE need this function here!!!
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Square ROI class
				
class squareROI(ROI):
	
	def __init__(self,embryo,name,Id,offset,sidelength,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.sidelength=sidelength
		self.offset=offset
		
	def setSideLength(self,s):
		self.sidelength=s
		return self.sidelength

	def getSideLength(self):
		return self.sidelength
	
	def setOffset(self,c):
		self.offset=c
		return self.offset
	
	def getOffset(self):
		return self.offset
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getSquareIdxImg(self.offset,self.sidelength,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getSquareIdxMesh(self.sidelength,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
		
		patch = ptc.Rectangle(self.offset,self.sidelength,self.sidelength,fill=False,linewidth=linewidth,color=color)
		ax.add_patch(patch)
		return ax
	
	def centerOffset(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return self.setOffset([self.embryo.dataResPx/2+0.5-self.sidelength/2,self.embryo.dataResPx/2+0.5-self.sidelength/2])
		else:
			return self.setOffset([self.embryo.dataResPx/2-self.sidelength/2,self.embryo.dataResPx/2-self.sidelength/2])
		
	def makeReducable(self,auto=False,debug=False):
		
		oldOffset=self.getOffset()
		self.centerOffset()
		
		if not auto:
			a=raw_input("Change offset of region "+ self.name + " from " + str(oldOffset) + ' to ' + str(self.getOffset()) + ' ? [Y/N]')
			
			if a=='N':
				self.setOffset(oldOffset)
				return False
			elif a=='Y':
				pass
			
			if not self.checkSymmetry():
				printWarning('Cannot make region '+self.name+' reducable.')
				self.setOffset(oldCenter)
				return False
			return True
		return False
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsideSquare(x,y,self.offset,self.sidelength)
	
	def computeXYExtend(self):
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelength]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelength]
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Square and slice ROI class

class squareSliceROI(squareROI,sliceROI):
	
	def __init__(self,embryo,name,Id,offset,sidelength,height,width,sliceBottom,color='b'):	
		
		squareROI.__init__(self,embryo,name,Id,offset,sidelength,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getSquareIdxImg(self.offset,self.sidelength,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getSquareIdxMesh(self.sidelength,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsideSquare(x,y,self.offset,self.sidelength)
	
	def computeXYExtend(self):
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelength]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelength]
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Rectangle ROI class

class rectangleROI(ROI):
	
	def __init__(self,embryo,name,Id,offset,sidelengthX,sidelengthY,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.sidelengthX=sidelengthX
		self.sidelengthY=sidelengthY
		self.offset=offset
		
	def setSideLengthX(self,s):
		self.sidelengthX=s
		return self.sidelengthX

	def getSideLengthX(self):
		return self.sidelengthX
	
	def setSideLengthY(self,s):
		self.sidelengthY=s
		return self.sidelengthY

	def getSideLengthY(self):
		return self.sidelengthY
	
	def setOffset(self,c):
		self.offset=c
		return self.offset
	
	def getOffset(self):
		return self.offset
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getRectangleIdxImg(self.offset,self.sidelengthX,self.sidelengthY,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getRectangleIdxMesh(self.sidelength,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
			
		patch = ptc.Rectangle(self.offset,self.sidelengthX,self.sidelengthY,fill=False,linewidth=linewidth,color=color)
		ax.add_patch(patch)
		return ax
	
	def centerOffset(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return self.setOffset([self.embryo.dataResPx/2+0.5-self.sidelengthX/2,self.embryo.dataResPx/2+0.5-self.sidelengthY/2])
		else:
			return self.setOffset([self.embryo.dataResPx/2-self.sidelengthX/2,self.embryo.dataResPx/2-self.sidelengthY/2])
		
	def makeReducable(self,atuo=False,debug=False):
		
		oldOffset=self.getOffset()
		self.centerOffset()
		
		if not auto:
			a=raw_input("Change offset of region "+ self.name + " from " + str(oldOffset) + ' to ' + str(self.getOffset()) + ' ? [Y/N]')
			
			if a=='N':
				self.setOffset(oldOffset)
				return False
			elif a=='Y':
				pass
			
			if not self.checkSymmetry():
				printWarning('Cannot make region '+self.name+' reducable.')
				self.setOffset(oldCenter)
				return False
			return True
		return False
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsideRectangle(x,y,self.offset,self.sidelengthX,self.sidelengthY)
	
	def computeXYExtend(self):
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelengthX]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelengthY]
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Rectangle and slice ROI class

class rectangleSliceROI(rectangleROI,sliceROI):
	
	def __init__(self,embryo,name,Id,offset,sidelength,height,width,sliceBottom,color='b'):
		
		rectangleROI.__init__(self,embryo,name,Id,offset,sidelengthX,sidelengthY,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getRectangleIdxImg(self.sidelengthX,self.sidelengthY,self.offset,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getRectangleIdxMesh(self.sidelengthX,self.sidelengthY,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx	
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsideRectangle(x,y,self.offset,self.sidelengthX,self.sidelengthY)
	
	def computeXYExtend(self):
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelengthX]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelengthY]
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Polygon ROI class
			
class polyROI(ROI):
	
	def __init__(self,embryo,name,Id,corners,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.corners=corners
	
	def setCorners(self,corners):
		self.corners=corners
		return corners
	
	def getCorners(self):
		return self.corners
	
	def addCorner(self,c,pos=-1):
		self.corners.insert(pos,c)
		return self.corners
	
	def appendCorner(self,c):
		self.corners.append(c)
		return self.corners
	
	def removeCorner(self,pos):
		self.corners.pop(pos)
		return self.corners
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getPolyIdxImg(self.corners,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getPolyIdxMesh(self.corners,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
		patch = ptc.Rectangle(self.corners,closed=True,fill=False,linewidth=linewidth,color=color)
		
		ax.add_patch(patch)
		return ax
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsidePoly(x,y,self.corners)
	
	def computeXYExtend(self):
		
		cornersNP=np.array(corners)
	
		xmax=cornersNP[:,0].max()
		xmin=cornersNP[:,0].min()
		ymax=cornersNP[:,1].max()
		ymin=cornersNP[:,1].min()
		
		self.xExtend=[xmin,xmax]
		self.yExtend=[ymin,ymax]
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Polygon and slice ROI class
			
class polySliceROI(polyROI,sliceROI):
	
	def __init__(self,embryo,name,Id,corners,height,width,sliceBottom,color='b'):	
		
		polyROI.__init__(self,embryo,name,Id,corners,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
		
		self.corners=corners
	
	def computeImgIdx(self,debug=False):
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getPolyIdxImg(self.corners,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		self.meshIdx=pyfrp_idx_module.getPolyIdxMesh(self.corners,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx	
	
	def checkXYInside(self,x,y):
		return pyfrp_idx_module.checkInsidePoly(x,y,self.corners)
	
	def computeXYExtend(self):
		
		cornersNP=np.array(corners)
	
		xmax=cornersNP[:,0].max()
		xmin=cornersNP[:,0].min()
		ymax=cornersNP[:,1].max()
		ymin=cornersNP[:,1].min()
		
		self.xExtend=[xmin,xmax]
		self.yExtend=[ymin,ymax]
		
		return self.xExtend, self.yExtend
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Polygon and slice ROI class
	
class customROI(ROI):
	
	def __init__(self,embryo,name,Id,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.ROIsIncluded=[]
		self.procedures=[]
		
	def addROI(self,r,p):
		if r not in self.ROIsIncluded:
			self.ROIsIncluded.append(r)
			self.procedures.append(p)
		return self.ROIsIncluded
	
	def removeROI(self,r):
		if r in self.ROIsIncluded:
			i=self.ROIsIncluded.index(r)
			self.ROIsIncluded.remove(r)
			self.procedures.pop(i)
			
		return self.ROIsIncluded
	
	def mergeROIs(self,r):
		
		if len(self.ROIsIncluded)==0:
			self.copyIdxs(r)
		else:
			self.computeImgMask()
			self.imgMask=self.imgMask*r.computeImgMask()
			
			self.imgIdxX,self.imgIdxY=pyfrp_idx_module.mask2ind(self.imgMask,self.embryo.dataResPx)
			
			self.meshIdx=pyfrp_misc_module.matchVals(self.meshIdx,r.meshIdx)
				
			self.addROI(r,1)
			
			self.extImgIdxX,self.extImgIdxY = pyfrp_idx_module.getCommonExtendedPixels(self.ROIsIncluded,self.embryo.dataResPx,procedures=self.procedures,debug=False)
			
		return self.ROIsIncluded
		
	def substractROIs(self,r):
		
		if len(self.ROIsIncluded)==0:
			self.copyIdxs(r)
		else:
			self.computeImgMask()
			self.imgMask=self.imgMask*(1-r.computeImgMask())
			self.imgIdxX,self.imgIdxY=pyfrp_idx_module.mask2ind(self.imgMask,self.embryo.dataResPx)
			
			self.meshIdx=pyfrp_misc_module.complValsSimple(self.meshIdx,r.meshIdx)
				
			self.addROI(r,-1)
			
			self.extImgIdxX,self.extImgIdxY = pyfrp_idx_module.getCommonExtendedPixels(self.ROIsIncluded,self.embryo.dataResPx,procedures=self.procedures,debug=False)
			
		return self.ROIsIncluded
	
	def getROIsIncluded(self):
		return self.ROIsIncluded
		
	def setROIsIncluded(self,l):
		self.ROIsIncluded=l
		return self.ROIsIncluded
	
	def updateIdxs(self):
			
		self.emptyIdxs()

		for i in range(len(self.ROIsIncluded)):
			if i==0:
				self.copyIdxs(self.ROIsIncluded[i])
			else:
				if self.procedures[i]==1:
					self.mergeROIs(self.ROIsIncluded[i])
				elif self.procedures[i]==-1:
					self.substractROIs(self.ROIsIncluded[i])
				else:
					printWarning("Unknown Procedure" + str(self.procedures[i]) + " in Custom ROI " + self.name +". Not going to do anything.")
					
		return self.getAllIdxs()
		
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		if color==None:
			color=self.color
		
		if color=='each':
			color=None
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
			
		for r in self.ROIsIncluded:
			if hasattr(r,'showBoundary'):
				r.showBoundary(color=color,ax=ax,linewidth=linewidth)
		
		return ax
	
	def checkXYInside(self,x,y):
		b=True
		for i,r in enumerate(self.ROIsIncluded):
			if self.procedures[i]==1:
				b=b and r.checkXYInside(x,y)
			elif self.procedures[i]==-1:
				b=b and not r.checkXYInside(x,y)
		return b
			
	def computeXYExtend(self):
		self.xExtend,self.yExtend=pyfrp_idx_module.getCommonXYExtend(self.ROIsIncluded)
		return self.xExtend,self.yExtend
	
	def roiIncluded(self,r):
		return r in self.ROIsIncluded
	
	