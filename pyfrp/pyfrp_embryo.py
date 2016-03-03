#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
import pyfrp_misc_module 
import pyfrp_img_module 
import pyfrp_fit_module 
import pyfrp_stats_module 
import pyfrp_IO_module
import pyfrp_plot_module
from pyfrp_term_module import *

#PyFRAP Objects
import pyfrp_geometry
import pyfrp_simulation
import pyfrp_analysis
import pyfrp_ROI
import pyfrp_fit


#matplotlib
import matplotlib.pyplot as plt

#Time 
import time

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Embryo object

class embryo:
	
	#Creates new embryo object
	def __init__(self,name):
		
		#Dataset name
		self.name = name
		
		#Data Image specifics
		self.dataEnc='uint16'
		self.dataFT='tif'
		self.dataResPx=512.
		self.dataResMu=322.34
		
		#DataSet files
		self.fileList=[]
		self.fnDatafolder=""
		
		#Data time specifics
		self.frameInterval=10
		self.nFrames=300
		self.tStart=0
		self.tEnd=self.tStart+self.frameInterval*(self.nFrames-1)
		self.tvecData=np.linspace(self.tStart,self.tEnd,self.nFrames)
		
		#Imaging specifics (In mu)
		self.sliceDepthMu=30.
		self.sliceWidthMu=2.5
		
		#Zoom Imaging (need to call this different)
		self.convFact=self.computeConvFact(updateDim=False)
		
		#Imaging specifics (In pixels)
		self.sliceWidthPx=self.sliceWidthMu/self.convFact
		self.sliceDepthPx=self.sliceDepthMu/self.convFact
		self.sliceHeightPx=-self.sliceDepthPx
		self.sliceBottom=0
		
		#Bleaching specifics (in Mu)
		self.sideLengthBleachedMu=2*79.54
		
		#Bleaching specifics (in pxs) (can take out)
		self.sideLengthBleachedPx=self.sideLengthBleachedMu/self.convFact
		self.offsetBleachedPx=[self.dataResPx/2-self.sideLengthBleachedPx/2, self.dataResPx/2-self.sideLengthBleachedPx/2]
		
		#List of ROIs
		self.ROIs=[]
		
		#Master ROI
		self.masterROIIdx=None
		
		#Geometry
		self.geometry=None
		
		#Simulation
		self.simulation=None
		
		#Analysis
		self.analysis=None
		
		#Fitting 
		self.fits=[]

	def addFit(self,fit):
		self.fits.append(fit)
		return self.fits
	
	def newFit(self,name):
		fit=pyfrp_fit.fit(self,name)
		self.addFit(fit)
		return fit
		
	def deleteFit(self,i):
		self.fits.pop(i)
		return self.fits
	
		
	def save(self,fn=None):
		
		if fn==None:
			fn=self.name+".emb"
		
		pyfrp_IO_module.saveToPickle(self,fn=fn)	
		print "Saved "+  self.name+ " to " + fn
		return fn
		
	def copy(self):
		copiedEmbryo=cpy.deepcopy(self)
		return copiedEmbryo
	
	def updateVersion(self):
		embtemp=embryo("temp")
		pyfrp_misc_module.updateObj(embtemp,self)
		return self
	
	def computeConvFact(self,updateDim=True):
		self.convFact = self.dataResMu/self.dataResPx
		if updateDim:
			self.updatePxDimensions()
		return self.convFact
	
	def updatePxDimensions(self):
		
		self.sliceWidthPx=self.sliceWidthMu/self.convFact
		self.sliceDepthPx=self.sliceDepthMu/self.convFact
		self.sliceHeightPx=-self.sliceDepthPx
			
		self.updateBleachedRegion()

		return self.sliceWidthPx,self.sliceDepthPx, self.sliceHeightPx, self.sideLengthBleachedPx, self.offsetBleachedPx
	
	def setName(self,n):
		self.name=n
		return self.name
	
	def getName(self):
		return self.name
	
	def setDataResMu(self,res):
		self.dataResMu=res 
		self.computeConvFact()
		return self.dataResMu
	
	def setDataResPx(self,res):
		self.dataResPx=res 
		self.computeConvFact()
		return self.dataResPx
	
	def setDataFT(self,f):
		self.dataFT=f 
		return self.dataFT
	
	def setDataEnc(self,e):
		self.dataEnc=e 
		return self.dataEnc
	
	def getDataFT(self):
		return self.dataFT
	
	def getDataEnc(self):
		return self.dataEnc
	
	def setFrameInterval(self,dt):
		self.frameInterval=dt
		self.updateTimeDimensions()
		return self.frameInterval
	
	def setTStart(self,t):
		self.tStart=t
		self.updateTimeDimensions()
		return self.tStart
	
	def setTEnd(self,t):
		self.tEnd=t
		self.updateTimeDimensions()
		return self.tEnd
	
	def updateNFrames(self):
		self.nFrames=len(self.fileList)
		return self.nFrames
	
	def getNFrames(self):
		return self.nFrames
	
	def setNFrames(self,n):
		self.nFrames=n
		return self.nFrames
	
	def getDataResPx(self):
		return self.dataResPx
	
	def getDataResMu(self):
		return self.dataResMu
	
	def getFrameInterval(self):
		return self.frameInterval
	
	def getTStart(self):
		return self.tStart
	
	def getTEnd(self):
		return self.tEnd
	
	def getTvecData(self):
		return self.tvecData
	
	def updateTimeDimensions(self):
		self.tEnd=self.tStart+self.frameInterval*(self.nFrames-1)
		self.tvecData=np.linspace(self.tStart,self.tEnd,self.nFrames)
		if self.simulation!=None:
			self.simulation.toDefaultTvec()
		return self.tvecData
	
	def setSliceDepthMu(self,d):
		self.sliceDepthMu=d
		self.updatePxDimensions()
		return self.sliceDepthMu
	
	def setMasterROIIdx(self,r):
		self.masterROI=r
		return self.masterROI
	
	def setSideLengthBleachedMu(self,s):
		self.sideLengthBleachedMu=s
		self.updateBleachedRegion()
		return self.sideLengthBleachedMu
		
	def updateBleachedRegion(self):
		self.sideLengthBleachedPx=self.sideLengthBleachedMu/self.convFact
		self.offsetBleachedPx=[self.dataResPx/2-self.sideLengthBleachedPx/2, self.dataResPx/2-self.sideLengthBleachedPx/2]
		return self.sideLengthBleachedPx, self.offsetBleachedPx
		
	def getMasterROIIdx(self):
		return self.masterROIIdx
	
	def newROI(self,name,Id,zmin=-np.inf,zmax=np.inf,color='b',asMaster=False):
		roi=pyfrp_ROI.ROI(self,name,Id,zmin=-np.inf,zmax=np.inf,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
		
	def newRadialROI(self,name,Id,center,radius,color='b',asMaster=False):
		roi=pyfrp_ROI.radialROI(self,name,Id,center,radius,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
		
	def newSliceROI(self,name,Id,height,width,sliceBottom,color='b',asMaster=False):
		roi=pyfrp_ROI.sliceROI(self,name,Id,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newRadialSliceROI(self,name,Id,center,radius,height,width,sliceBottom,color='b',asMaster=False):
		roi=pyfrp_ROI.radialSliceROI(self,name,Id,center,radius,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newSquareROI(self,name,Id,offset,sidelength,color='b',asMaster=False):
		roi=pyfrp_ROI.squareROI(self,name,Id,offset,sidelength)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newSquareSliceROI(self,name,Id,offset,sidelength,height,width,sliceBottom,color='b',asMaster=False):
		roi=pyfrp_ROI.squareSliceROI(self,name,Id,offset,sidelength,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newRectangleROI(self,name,Id,offset,sidelengthX,sidelengthY,color='b',asMaster=False):
		roi=pyfrp_ROI.rectangleROI(self,name,Id,offset,sidelengthX,sidelengthY,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newRectangleSliceROI(self,name,Id,offset,sidelength,height,width,sliceBottom,color='b',asMaster=False):
		roi=pyfrp_ROI.rectangleSliceROI(self,name,Id,offset,sidelength,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newPolyROI(self,name,Id,corners,color='b',asMaster=False):
		roi=pyfrp_ROI.polyROI(self,name,Id,corners,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newPolySliceROI(self,name,Id,corners,height,width,sliceBottom,color='b',asMaster=False):
		roi=pyfrp_ROI.polySliceROI(self,name,Id,corners,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newCustomROI(self,name,Id,color='b',asMaster=False):
		roi=pyfrp_ROI.customROI(self,name,Id,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def getFreeROIId(self):
		IdList=pyfrp_misc_module.objAttrToList(self.ROIs,'Id')
		notFound=True
		i=0
		while notFound:
			if i not in IdList:
				return i
			i=i+1
			
	def removeROI(self,i):
		self.ROIs.pop(i)
		if i==self.masterROIIdx:
			printWarning('You are deleting the masterROI. You need to select a new one before you can analyze or simulate.')
		return self.ROIs
	
	def listROIs(self):
		for r in self.ROIs:
			print r.name , type(r)
		return True
	
	def genDefaultROIs(self,center,radius,rimFactor=0.66):
		
		allsl=self.newSliceROI("All",0,-np.inf,np.inf,True,color=(0.5,0.3,0.4))
		
		allsqu=self.newSquareROI('All Square',1,self.offsetBleachedPx,self.sideLengthBleachedPx,color=(0.1,0.7,0.5))
		
		allout=self.newCustomROI("All Out",2,color=(0.3,0.4,0.5))
		allout.addROI(allsl,1)
		allout.addROI(allsqu,-1)
		
		sl=self.newRadialSliceROI("Slice",3,center,radius,self.sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='g',asMaster=True)
		
		sl2=self.newRadialSliceROI("Slice rim",6,center,rimFactor*radius,self.sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='y')
		 
		rim=self.newCustomROI("Rim",7,color='m')
		rim.addROI(sl,1)
		rim.addROI(sl2,-1)
		rim.setUseForRim(True)
		
		squ=self.newSquareSliceROI("Bleached Square",4,self.offsetBleachedPx,self.sideLengthBleachedPx,self.sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='b')
		
		out=self.newCustomROI("Out",5,color='r')
		out.addROI(sl,1)
		out.addROI(squ,-1)
		
		return self.ROIs
	
	def getROIs(self):
		return self.ROIs
	
	def getROIByName(self,name):
		for r in self.ROIs:
			if r.name==name:
				return r
		
		printWarning('Cannot find ROI with name '+ name + '. Will return None. This can lead to further problems.')
		
		return None

	def getROIById(self,Id):
		for r in self.ROIs:
			if r.Id==Id:
				return r
			
		
		printWarning('Cannot find ROI with Id '+ Id + '. Will return None. This can lead to further problems.')
		
		return None
	
	def updateFileList(self):
		self.fileList=pyfrp_misc_module.getSortedFileList(self.fnDatafolder,self.dataFT)
		if len(self.fileList)<0:
			printWarning("There are no files of type " + self.dataFT + " in " + self.fnDatafolder +" . This can lead to problems.")
		return self.fileList
	
	def setFileList(self,l):
		self.fileList=l
		return self.fileList
	
	def getFileList(self):
		return self.fileList
	
	def setDataFolder(self,fn):
		self.fnDatafolder=fn
		return self.fnDatafolder
	
	def getDataFolder(self):
		return self.fnDatafolder
	
	def getGeometry(self):
		return self.geometry
	
	def setGeometry2ZebraFishDomeStage(self,center,imagingRadius,radiusScale=1.1):
		self.geometry=pyfrp_geometry.zebrafishDomeStage(self,center,imagingRadius,radiusScale=radiusScale)
		return self.geometry
	
	def setGeometry2Cylinder(self,center,radius,height):
		self.geometry=pyfrp_geometry.cylinder(self,center,radius,height)
		return self.geometry
	
	def setGeometry2Cone(self,center,upperRadius,lowerRadius,height):
		self.geometry=pyfrp_geometry.cone(self,center,upperRadius,lowerRadius,height)
		return self.geometry
	
	def setGeometry2Ball(self,center,imagingRadius):
		self.geometry=pyfrp_geometry.xenopusBall(self,center,imagingRadius)
		return self.geometry
	
	def setGeometry2ZebraFishDomeStageQuad(self,center,imagingRadius,radiusScale=1.1):
		self.geometry=pyfrp_geometry.zebrafishDomeStageQuad(self,center,imagingRadius,radiusScale=radiusScale)
		return self.geometry
	
	def setGeometry2CylinderQuad(self,center,radius,height):
		self.geometry=pyfrp_geometry.cylinderQuad(self,center,radius,height)
		return self.geometry
	
	def setGeometry2BallQuad(self,center,imagingRadius):
		self.geometry=pyfrp_geometry.xenopusBallQuad(self,center,imagingRadius)
		return self.geometry
	
	def setGeometry2Custom(self,center,fnGeo=""):
		self.geometry=pyfrp_geometry.custom(self,center,fnGeo)
		return self.geometry
	
	def newSimulation(self):
		self.simulation=pyfrp_simulation.simulation(self)
		return self.simulation
	
	def newAnalysis(self):
		self.analysis=pyfrp_analysis.analysis(self)
		return self.analysis
	
	def computeROIIdxs(self,signal=None,debug=True):
			
		for i,r in enumerate(self.ROIs):
			startInit=time.clock()
			r.computeIdxs()
			
			if debug:
				print r.name, time.clock()-startInit
			if signal:
				signal.emit(int(100.*i)/float(len(self.ROIs)))
			
		return self.ROIs	
	
	def geometry2Quad(self):
	
		if self.geometry==None:
			printError("No geometry selected yet.")
			return
		
		if 'Quad' in self.geometry.typ:
			printError("Geometry already set to quad.")
			return
		
		if self.geometry.typ=='zebrafishDomeStage':
			self.setGeometry2BallQuad(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		elif self.geometry.typ=='cylinder':
			self.setGeometry2CylinderQuad(self,self.geometry.getCenter(),self.geomtry.getRadius())
		elif self.geometry.typ=='xenopusBall':
			self.setGeometry2BallQuad(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		
		return self.geometry	
			
	
	def geometry2Full(self):
		
		if self.geometry==None:
			printError("No geometry selected yet.")
			return
		
		if 'Quad' not in self.geometry.typ:
			printError("Geometry already set to full.")
			return
		
		if self.geometry.typ=='zebrafishDomeStageQuad':
			self.setGeometry2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		elif self.geometry.typ=='cylinderQuad':
			self.setGeometry2Cylinder(self,self.geometry.getCenter(),self.geomtry.getRadius())
		elif self.geometry.typ=='xenopusBallQuad':
			self.setGeometry2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		
		return self.geometry	
	
	def setEmbryo2Quad(self):	
		self.geometry2Quad()
		self.ROIs2Quad()
		return self.geometry,self.ROIs
	
	def setEmbryo2Full(self):	
		self.geometry2Full()
		self.ROIs2Full()
		return self.geometry,self.ROIs
	
	def ROIs2Quad(self):
		for r in self.ROIs:
			r.idxs2Quad()
		return self.ROIs
	
	def ROIs2Full(self):
		for r in self.ROIs:
			r.idxs2Full()
		return self.ROIs
	
	def showAllROIBoundaries(self,ax=None):
		for r in self.ROIs:
			if hasattr(r,'showBoundary'):
				ax=r.showBoundary(ax=ax)
		return ax
	
	def loadDataImg(self,idx):
		return pyfrp_img_module.loadImg(self.fnDatafolder+self.fileList[idx],self.dataEnc)
	
	def showDataImg(self,ax=None,idx=0):
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],sup=["Embryo" + self.name],titles=["DataImg "+str(idx)],tight=False)
			ax=axes[0]
			
		img=self.loadDataImg(idx)	
		ax.imshow(img)
		
		return ax
		
	def showAllROIIdxs(self,axes=None):
		
		if axes==None:
			proj=len(self.ROIs)*['3d']+len(self.ROIs)*[None]+len(self.ROIs)*[None]
			fig,axes = pyfrp_plot_module.makeSubplot([3,len(self.ROIs)],sup=["Embryo" + self.name + " ROI Indices"],tight=True,proj=proj)
			
		for i,r in enumerate(self.ROIs):
			print r.name
			currAxes=[axes[0+i],axes[len(self.ROIs)+i],axes[2*len(self.ROIs)+i]]
			r.showIdxs(axes=currAxes)
			
		return axes
			
	
	def plotAllData(self,ax=None,legend=True):
		for r in self.ROIs:
			ax=r.plotData(ax=ax)
		
		#if legend:
			#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		return ax	
	
	def plotAllSim(self,ax=None,legend=True):
		for r in self.ROIs:
			ax=r.plotSim(ax=ax)
		
		if legend:
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		return ax	
	
	def plotAllDataPinned(self,ax=None,legend=True):
		for r in self.ROIs:
			ax=r.plotDataPinned(ax=ax)
		
		if legend:
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		return ax	
	
	def plotAllSimPinned(self,ax=None,legend=True):
		for r in self.ROIs:
			ax=r.plotSimPinned(ax=ax)
		
		if legend:
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		return ax	
	
	
	def checkQuadReducable(self,tryFix=False,auto=False,debug=False):
		reducable=True
		for r in self.ROIs:
			if not r.checkSymmetry(debug=debug):
				if tryFix:
					if hasattr(r,'makeReducable'):
						reducable=r.makeReducable(auto=auto,debug=debug)
				else:		
					printWarning('Cannot reduce region '+ self.name+' to quadrant. Indices are not symmetric.')
					reducable=False
				
		return reducable		
	
	def makeQuadReducable(self,auto=False,debug=False):
		reducable=self.checkQuadReducable(tryFix=True,auto=auto,debug=debug)

		if reducable:
			self.geometry.centerInImg()
		return reducable
	
	def getMasterROI(self):
		return self.ROIs[self.masterROIIdx]
	
	def getOptimalAllROI(self,name='All',makeNew=False):
		self.geometry.setAllROI(name=name,makeNew=makeNew)
	
	def computePinVals(self,useMin=True,useMax=True,bkgdVal=None,debug=False):
		
		if bkgdVal==None:
			bkgdVal=self.computeBkgd(useMin=useMin,debug=debug)
		normVal=self.computeNorm(bkgdVal,useMax=useMax,debug=debug)
		
		return bkgdVal,normVal
	
	def computeBkgd(self,useMin=False,fromTS='both',debug=False):
		bkgds=[]
		for r in self.ROIs:
			if fromTS in ['data','both']:
				bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.dataVec,useMin=useMin,useMax=False,debug=debug)
				bkgds.append(bkgdTemp)
			if fromTS in ['sim','both']:
				bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.simVec,useMin=useMin,useMax=False,debug=debug)
				bkgds.append(bkgdTemp)
			
			print r.name,bkgds[-2:]
				
		return min(bkgds)	
	
	def computeNorm(self,bkgdVal,useMax=True,fromTS='both',debug=False):
		norms=[]
		for r in self.ROIs:
			if fromTS in ['data','both']:
				bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.dataVec,bkgdVal=bkgdVal,useMin=False,useMax=useMax,debug=debug)
				norms.append(norm)
			if fromTS in ['sim','both']:
				bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.simVec,bkgdVal=bkgdVal,useMin=False,useMax=useMax,debug=debug)
				norms.append(norm)
		return max(norms)	
			
	def computeIdealFRAPPinVals(self,bkgdName='Bleached Square',normName='Slice',debug=False,useMin=False,useMax=False,sepSim=True,switchThresh=0.95):
	
		bkgdROI=self.getROIByName(bkgdName)
		normROI=self.getROIByName(normName)
		
		if bkgdROI==None or normROI==None:
			return None,None,None,None
		
		bkgdValSim=None
		normValSim=None
		
		bkgds=[]
		norms=[]
		
		bkgdBdata,normBdata=pyfrp_fit_module.computePinVals(bkgdROI.dataVec,useMin=useMin,useMax=useMax,debug=debug)
		bkgdBsim,normBsim=pyfrp_fit_module.computePinVals(bkgdROI.simVec,useMin=useMin,useMax=useMax,debug=debug)
		
		bkgdVal=min([bkgdBdata,bkgdBsim])
		
		if sepSim:
			bkgdValSim=bkgdBsim
			bkgdVal=bkgdBdata
		
		"""Need to check if we need to switch ROI. This is important if recovery curves are very slow and bleached region
		does not reach full recovery. If this happens, we divide by a number that is smaller than 1 and hence boost all 
		timeseries way above one instead of limiting it below one."""
		
		if normBdata<switchThresh or normBsim<switchThresh:
			bkgdNdata,normNdata=pyfrp_fit_module.computePinVals(normROI.dataVec,useMin=useMin,useMax=useMax,bkgdVal=bkgdVal,debug=debug)
			bkgdNsim,normNsim=pyfrp_fit_module.computePinVals(normROI.simVec,useMin=useMin,useMax=useMax,bkgdVal=bkgdVal,debug=debug)
			
			normVal=max([normNsim,normNdata])
			if sepSim:
				normValSim=normNsim
				normVal=normNdata
			
			if debug:
				printNote('Switched to region ' +  normName + ' for computation of normalization value.' )
		else:	
			normVal=max([normBdata,normBsim])
			if sepSim:
				normValSim=normBsim
				normVal=normBdata
		
		return bkgdVal, normVal, bkgdValSim, normValSim
		
	def pinAllROIs(self,bkgdVal=None,normVal=None,bkgdValSim=None,normValSim=None,useMin=False,useMax=False,debug=False):
		
		bkgdValTemp,normValTemp = self.computePinVals(useMin=useMin,useMax=useMax,debug=debug)
		
		bkgdVal = pyfrp_misc_module.assignIfVal(bkgdVal,bkgdValTemp,None)
		normVal = pyfrp_misc_module.assignIfVal(normVal,normValTemp,None)
		
		for r in self.ROIs:
			r.pinAllTS(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=debug)
		return self.ROIs	
			
	def printAllAttr(self):
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Dont print out large arrays
			if isinstance(vars(self)[str(item)],(int,float,str)):
				print item, " = ", vars(self)[str(item)]
			elif isinstance(vars(self)[str(item)],(list,np.ndarray)) and max(np.shape(vars(self)[str(item)]))<5:
				print item, " = ", vars(self)[str(item)]		
	
	def isAnalyzed(self):
		b=True
		for r in self.ROIs:
			b=b and r.isAnalyzed()
		return b

	def isSimulated(self):
		b=True
		for r in self.ROIs:
			b=b and r.isSimulated()
		return b
	
	def isFitted(self):
		b=True
		for fit in self.fits:
			b=b and fit.isFitted()
		return b	
			
	def getFitByName(self,name):
		for fit in self.fits:
			if fit.name==name:
				return fit
			
		printWarning('Cannot find fit with name '+ name + '. Will return None. This can lead to further problems.')
		
		return None
	
	def getInterpolationError(self):
		for r in self.ROIs:
			print r.name, r.getInterpolationError()
		
	
	#def grabDataDetails(self):
	###NOTE make a function that automatically grabs filetype and what not
	
	