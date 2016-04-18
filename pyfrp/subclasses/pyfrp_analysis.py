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
from pyfrp.modules import pyfrp_img_module 
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules.pyfrp_term_module import *

#matplotlib
import matplotlib.pyplot as plt

#===========================================================================================================================================================================
#Module Classes
#===========================================================================================================================================================================


class analysis:
	
	"""PyFRAP analysis class storing information about analysis options 
	and some analysis results.
	
	Analysis options are:
	
		* ``gaussian``: Apply gaussian filter to images. Default kernel size is ``gaussianSigma=2``.
		* ``median``: Apply gaussian filter to images. Default kernel size is ``medianRadius=5``.
		* ``flatten``: Apply flattening mask.
		* ``norm``: Norm  by pre image.
		* ``bkgd``: Substract background.
		* ``quad``: Perform reduction to first quadrant by flipping.
		* ``flipBeforeProcess``: Flip into quadrant before other processing options are applied.
	
	Analysis options are stored in ``process`` dictionary. If analysis finds option in ``process.keys``, it will
	perform option. Analysis options can be turned on/off using the respective functions, such as 
	
		* :py:func:pyfrp.subclasses.pyfrp_analysis.medianOn
		* :py:func:pyfrp.subclasses.pyfrp_analysis.flattenOn
		* etc.
	
	Processing parameters are stored in ``process.values``.
	
	The default processing options are ``process={}``, meaning that no image modification is applied before 
	concentration readout.
	
	.. warning: Quadrant reduction is still experimental.	
		
	Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): PyFRAP embryo instance.
	
	"""
	
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
		
		"""Runs analysis by passing analysis object to :py:func:`pyfrp.modules.pyfrp_img_module.analyzeDataset`.
		
		Will first check if ROI indices are computed for all ROIs and if necessary compute them before starting 
		data analysis.
		
		Keyword Args:
			signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
			embCount (int): Counter of counter process if multiple datasets are analyzed. 
			debug (bool): Print final debugging messages and show debugging plots.
			debugAll (bool): Print debugging messages and show debugging plots of each step.
			showProgress (bool): Print out progress.
		
		Returns:
			pyfrp.subclasses.pyfrp_analysis.analysis: Updated analysis instance.
		
		"""
		
		if 'norm' in self.process and 'flatten' in  self.process:
			printWarning("Both norm and flatten have been selected for data analysis. This is not advisable.")
		
		if not self.embryo.checkROIIdxs()[0]:
			self.embryo.computeROIIdxs()
			
		self=pyfrp_img_module.analyzeDataset(self,signal=signal,embCount=embCount,debug=debug,debugAll=debugAll,showProgress=showProgress)
		return self
	
	def setGaussianSigma(self,s):
		
		"""Sets size of gaussian kernel and updates its value
		in ``process`` dictionary if gaussian filter is turned on.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian_filter.
		
		Args:
			s (float): New sigma.
			
		"""
		
		self.gaussianSigma=s
		self.updateProcess()
		return self.gaussianSigma
	
	def getGaussianSigma(self):
		
		"""Returns size of gaussian kernel.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian_filter.
		
		Returns:
			float: Gaussian sigma.
			
		"""
		
		return self.gaussianSigma
	
	def setMedianRadius(self,s):
		
		"""Sets size of median kernel and updates its value
		in ``process`` dictionary if median filter is turned on.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.median and 
		http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.ndimage.filters.median_filter.html.
		
		Args:
			s (float): New radius.
			
		"""
		
		self.medianRadius=s
		self.updateProcess()
		return self.medianRadius
	
	def getMedianRadius(self):
		
		"""Returns size of median kernel.
		
		See also http://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.median and 
		http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.ndimage.filters.median_filter.html.
		
		Returns:
			float: New radius.
			
		"""
		
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
		
		self.flatteningMask=pyfrp_img_module.computeFlatMask(meanImg,self.dataOffset)
	
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
		
		fileList=pyfrp_misc_module.getSortedFileList(self.fnPreimage,self.embryo.dataFT)
		fileList=fileList[:self.nPre]
		meanImg=pyfrp_img_module.computeMeanImg(self.fnPreimage,fileList,self.embryo.dataEnc)
		
		if applyProcess:
		
			processDic=dict(self.process)
			processDic=self.removeProcessStep(processDic,'norm')
			
			imgPre=pyfrp_img_module.processImg(meanImg,processDic,flatteningMask,bkgdMask,None,dataOffset=self.dataOffset)
		
		self.preMask=imgPre
		
		return imgPre
 
	def getOptimalOffset(self,debug=False):
		
		self.dataOffset=pyfrp_img_module.findMinOffset(self.embryo.fnDatafolder,self.embryo.fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
		
		if self.fnPreimage!=None:
			fileList=pyfrp_misc_module.getSortedFileList(self.fnPreimage,self.embryo.dataFT)[:self.nPre]
			self.dataOffset=pyfrp_img_module.findMinOffset(self.fnPreimage,fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
		
		if self.fnFlatten!=None:
			fileList=pyfrp_misc_module.getSortedFileList(self.fnFlatten,self.embryo.dataFT)[:self.nFlatten]
			self.dataOffset=pyfrp_img_module.findMinOffset(self.fnFlatten,fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
			
		if self.fnBkgd!=None:
			fileList=pyfrp_misc_module.getSortedFileList(self.fnBkgd,self.embryo.dataFT)[:self.nBkgd]
			self.dataOffset=pyfrp_img_module.findMinOffset(self.fnBkgd,fileList,self.embryo.dataEnc,oldOffset=self.dataOffset,defaultAdd=1.,debug=debug)
			
		return self.dataOffset
		
	###NOTE: This function can be removed at publication date
	def testImport(self):
		printWarning("Testing module")
		print pyfrp_misc_module.remRepeatsList([3,3,4,56,1,3])
		
		