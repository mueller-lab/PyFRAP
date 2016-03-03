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
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
from numpy import *

#JSON Compression
import pickle

#PyFRAP Modules
import pyfrp_img_module 
import pyfrp_misc_module 
from pyfrp_term_module import *

#matplotlib
import matplotlib.pyplot as plt

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#analysis object

class analysis:
	
	#Creates new embryo object
	def __init__(self,embryo):
		
		self.embryo=embryo
		
		#Rim handling
		self.addRimImg=False
		self.concRim=None
		
		#Norming Data
		self.fnPreimage=''
		self.fnFlatten=''
		self.fnBkgd=''
		self.flatteningMask=None
		self.bkgdMask=None
		self.preMask=None
		self.nPre=2
		self.nFlatten=2
		self.nBkgd=2
		
		#Data offset
		self.dataOffset=1.
			
		#Processing options
		self.process={}
		self.gaussianSigma=2
		self.medianRadius=5
		
	def run(self,signal=None,embCount=None,debug=False,debugAll=False,showProgress=True):
		if 'norm' in self.process and 'flatten' in  self.process:
			printWarning("Both norm and flatten have been selected for data analysis. This is not advisable.")

		self=pyfrp_img_module.analyzeDataset(self,signal=signal,embCount=embCount,debug=debug,debugAll=debugAll,showProgress=showProgress)
		return self
	
	def setGaussianSigma(self,s):
		self.gaussianSigma=s
		self.updateProcess()
		return self.gaussianSigma
	
	def getGaussianSigma(self):
		return self.gaussianSigma
	
	def setMedianRadius(self,s):
		self.medianRadius=s
		self.updateProcess()
		return self.medianRadius
	
	def getMedianRadius(self):
		return self.medianRadius
	
	def setGaussian(self,b):
		return self.parm2Process(b,'gaussian',self.gaussianSigma)
		
	def setNorm(self,b):
		return self.parm2Process(b,'norm',self.fnPreimage)
	
	def setFlipBeforeProcess(self,b):
		return self.parm2Process(b,'flipBeforeProcess',True)
		
	def setMedian(self,b):
		return self.parm2Process(b,'median',self.medianRadius)
	
	def setQuad(self,b):
		return self.parm2Process(b,'quad',True)
	
	def setFlatten(self,b):
		return self.parm2Process(b,'flatten',self.fnFlatten)
	
	def setBkgd(self,b):
		return self.parm2Process(b,'bkgd',self.fnBkgd)
	
	def setFnFlatten(self,fn):
		self.fnFlatten=fn
		self.updateProcess()
		return self.fnFlatten
	
	def getFnFlatten(self):
		return self.fnFlatten
	
	def setFnBkgd(self,fn):
		self.fnBkgd=fn
		self.updateProcess()
		return self.fnBkgd
	
	def getFnBkgd(self):
		return self.fnBkgd
	
	def setFnPre(self,fn):
		self.fnPreimage=fn
		self.updateProcess()
		return self.fnPreimage
	
	def getFnPre(self):
		return self.fnPreimage
	
	def setDataOffset(self,s):
		self.dataOffset=s 
		return self.dataOffset
	
	def getDataOffset(self):
		return self.dataOffset
	
	def setAddRimImg(self,s):
		self.addRimImg=s 
		return self.addRimImg
	
	def getAddRimImg(self):
		return self.addRimImg
	
	def setConcRim(self,s):
		self.concRim=s 
		return self.concRim
	
	def getConcRim(self):
		return self.concRim
	
	def setProcess(self,s):
		self.process=s 
		return self.process
	
	def getProcess(self):
		return self.process
	
	def setNPre(self,n):
		self.nPre=n 
		return self.nPre
	
	def getNPre(self):
		return self.nPre
	
	def setNFlatten(self,n):
		self.nFlatten=n 
		return self.nFlatten
	
	def getNFlatten(self):
		return self.nFlatten
	
	def setNBkgd(self,n):
		self.nBkgd=n 
		return self.nBkgd
	
	def getNBkgd(self):
		return self.nBkgd
	
	def normOn(self):
		return 'norm' in self.process.keys()
	
	def bkgdOn(self):
		return 'bkgd' in self.process.keys()
	
	def flattenOn(self):
		return 'flatten' in self.process.keys()
	
	def gaussianOn(self):
		return 'gaussian' in self.process.keys()
	
	def medianOn(self):
		return 'median' in self.process.keys()
	
	def quadOn(self):
		return 'quad' in self.process.keys()
	
	def flipBeforeProcessOn(self):
		return 'flipBeforeProcess' in self.process.keys()
	
	def updateProcess(self):
		if 'flatten' in self.process.keys():
			self.process['flatten']=self.fnFlatten
		if 'norm' in self.process.keys():
			self.process['norm']=self.fnPreimage
		if 'bkgd' in self.process.keys():
			self.process['bkgd']=self.fnBkgd	
		if 'gaussian' in self.process.keys():
			self.process['gaussian']=self.gaussianSigma
		if 'median' in self.process.keys():
			self.process['median']=self.medianRadius	
			
		return
	
	def printProcess(self):
		printDict(self.process)
		return
	
	def parm2Process(self,b,key,val):
		if b:
			self.process[key]=val
		else:	
			try:
				self.process.pop(key)
			except KeyError:
				pass
		return self.process	
		
	def genDefaultProcess(self):
		self.setGaussian(False)
		self.setFlipBeforeProcess(True)
		self.setQuad(False)
		self.setNorm(False)
		self.setMedian(False)
		self.setFlatten(False)
		self.setBkgd(False)
		return self.process
	
	def removeProcessStep(self,dic,step):
		try:
			dic.pop(step)
		except KeyError:
			pass
		return dic
		
	def computeFlatteningMask(self,applyProcess=True):
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnFlatten,self.embryo.dataFT)
		fileList=fileList[:self.nFlatten]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnFlatten,fileList,self.embryo.dataEnc)
		
		if applyProcess:
		
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			processDic=self.removeProcessStep(processDic,'flatten')
			processDic=self.removeProcessStep(processDic,'bkgd')
			
			###NOTE: Also remove bkgd??? Having bkgd in there could lead to singularities with flattening norming? 
			meanImg=pyfrp_img_module.processImg(meanImg,processDic,None,None,None,dataOffset=self.dataOffset)
		
		self.flatteningMask=pyfrp_img_module.computeFlatMask(meanImg)
	
		return self.flatteningMask
	
	def computeBkgdMask(self,flatteningMask,applyProcess=True,applyFlatten=False):
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnBkgd,self.embryo.dataFT)
		fileList=fileList[:self.nBkgd]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnBkgd,fileList,self.embryo.dataEnc)
		
		if applyProcess:
			
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			processDic=self.removeProcessStep(processDic,'bkgd')
			if not applyFlatten:
				processDic=self.removeProcessStep(processDic,'flatten')
			
			meanImg=pyfrp_img_module.processImg(meanImg,processDic,flatteningMask,None,None,dataOffset=self.dataOffset)
		
		self.bkgdMask=meanImg
		
		return self.bkgdMask
	
	def computePreMask(self,flatteningMask,bkgdMask,applyProcess=True):
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnPre,self.embryo.dataFT)
		fileList=fileList[:self.nPre]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnPre,fileList,self.embryo.dataEnc)
		
		if applyProcess:
		
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			
			imgPre=pyfrp_img_module.processImg(meanImg,processDic,flatteningMask,bkgdMask,None,dataOffset=self.dataOffset)
		
		self.preMask=imgPre
		
		return imgPre