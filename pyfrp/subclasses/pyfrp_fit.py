#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Fit module for PyFRAP toolbox, including following fit objects:

#(1) fit

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module 
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_fit_module

from pyfrp.modules.pyfrp_term_module import *

#Time 
import time

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main fit class
		
class fit:
	
	#Create new fit object
	def __init__(self,embryo,name):
		
		#General Settings
		self.name=name
		self.embryo=embryo
		
		#Optimization algorithm settings
		self.optMeth="Constrained Nelder-Mead"
		self.maxfun=1000
		self.optTol=1e-10
		
		#Dataseries selection
		self.ROIsFitted=[]
		
		#What parameters to fit
		self.fitProd=False
		self.fitDegr=False
		
		#Equalization and pinning
		self.equOn=True
		self.fitPinned=True
		
		#Intial guess
		self.x0=[10,0,0]
		
		#Bounds
		self.LBProd=0.
		self.UBProd=100.
		self.LBDegr=0.
		self.UBDegr=100.
		self.LBD=0.01
		self.UBD=150.
		self.bounds=None
		
		#More settings
		self.kineticTimeScale=100000.
		self.bruteInit=False		
		
		#Cutting tvec option
		self.fitCutOffT=False
		self.cutOffT=150
		self.cutOffStepSim=self.embryo.simulation.stepsSim
		self.cutOffStepData=self.embryo.nFrames
		
		#Fit tracking
		self.saveTrack=0
		self.trackedParms=[]
		self.trackedFits=[]
		
		#Fitted Vectors
		self.fittedVecs=[]
		self.tvecFit=None
		self.dataVecsFitted=[]
		
		#Results
		self.SSD=10000000
		self.DOptMu=None
		self.DOptPx=None
		self.prodOpt=None
		self.degrOpt=None
		self.success=None
		self.iterations=None
		self.fcalls=None
		
		#Statistics
		self.Rsq=None
		self.MeanRsq=None
		
		#Equalization
		self.equFacts=[]
		
		#Empty result dataseries
		self.tvecFit=embryo.tvecData

	def addROI(self,r):
		if r not in self.ROIsFitted:
			self.ROIsFitted.append(r)
		return self.ROIsFitted
	
	def addROIByName(self,name):
		r=self.embryo.getROIByName(name)
		return self.addROI(r)
	
	def addROIById(self,Id):
		r=self.embryo.getROIById(Id)
		return self.addROI(r)
	
	def removeROI(self,r):
		if r in self.ROIsFitted:
			self.ROIsFitted.remove(r)
		return self.ROIsFitted
	
	def getX0(self):
		
		#Copying x0 into local variable to pass to solver, pop entries that are currently not needed
		if self.fitProd and self.fitDegr:
			x0=list(self.x0)		
		elif self.fitProd and  not self.fitDegr:	
			x0=list(self.x0)
			x0.pop(2)
		elif not self.fitProd and self.fitDegr:
			x0=list(self.x0)
			x0.pop(1)		
		elif not self.fitProd and not self.fitDegr:
			x0=[self.x0[0]]
		
		return x0
	
	def getBounds(self):
	
		if self.fitProd and self.fitDegr:
			bnds = ((self.LBD, self.UBD), (self.LBD, self.UBD),(self.LBD,self.UBD))
			ranges=(slice(self.LBD,self.UBD,1),slice(self.LBProd,self.UBProd,10),slice(self.LBDegr,self.UBDegr,10))		
		elif self.fitProd and  not self.fitDegr:	
			bnds = ((self.LBD, self.UBD), (self.LBProd, self.UBProd))
			ranges=(slice(self.LBD,self.UBD,1),slice(self.LBProd,self.UBProd,10))
		elif not self.fitProd and self.fitDegr:
			bnds = ((self.LBD, self.UBD), (self.LBDegr, self.UBDegr))
			ranges=(slice(self.LBD,self.UBD,1),slice(self.LBDegr,self.UBDegr,10))
		elif not self.fitProd and not self.fitDegr:
			bnds = ((self.LBD, self.UBD),)
			ranges=(1,self.UBD)
			
		if self.optMeth=='brute':
			self.bounds=ranges
		else:
			self.bounds=bnds
			
		return self.bounds
		
	def run(self,debug=False,gui=None,ax=None):
		self=pyfrp_fit_module.FRAPFitting(self,debug=debug,gui=gui,ax=ax)
		return self
	
	def assignOptParms(self,res):
		if self.fitProd and self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=res[1]/self.kineticTimeScale
			self.degrOpt=res[2]/self.kineticTimeScale
			
		elif self.fitProd and not self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=res[1]/self.kineticTimeScale
			self.degrOpt=self.x0[2]/self.kineticTimeScale
			
		elif not self.fitProd and self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=self.x0[1]/self.kineticTimeScale
			self.degrOpt=res[1]/self.kineticTimeScale
			
		elif not self.fitProd and not self.fitDegr:
			self.DOptPx=res[0]
			self.prodOpt=self.x0[1]/self.kineticTimeScale
			self.degrOpt=self.x0[2]/self.kineticTimeScale
		
		
		self.DOptMu=self.DOptPx*self.embryo.convFact**2
		
		return self.DOptPx, self.prodOpt, self.degrOpt, self.DOptMu
	
	def plotFit(self,ax=None):
		
		for r in self.ROIsFitted:
			ax=r.plotFit(self,ax=ax)
			
		return ax
	
	def printResults(self):
		
		printObjAttr('DOptMu',self)
		printObjAttr('DOptPx',self)
		printObjAttr('prodOpt',self)
		printObjAttr('degrOpt',self)
		printObjAttr('success',self)
		printObjAttr('Rsq',self)
		printObjAttr('MeanRsq',self)
		
		return True
	
	def setOptMeth(self,m):
		self.optMeth=m
		return self.optMeth
	
	def getOptMeth(self):
		return self.optMeth
		
	def isFitted(self):
		return self.DOptMu!=None
		
	def setEqu(self,b):
		self.equOn=b
		return self.equOn
	
	def setFitPinned(self,b):
		self.fitPinned=b
		return self.fitPinned
	
	def setFitProd(self,b):
		self.fitProd=b
		return self.fitProd
	
	def setFitDegr(self,b):
		self.fitDegr=b
		return self.fitDegr
	
	def setSaveTrack(self,b):
		self.saveTrack=b
		return self.saveTrack
	
	def setFitCutOffT(self,b):
		printWarning("CutOffT Option is currently VERY experimental. Fitting might crash.")
		self.fitCutOffT=b
		return self.fitCutOffT
	
	def setCutOffT(self,t):
		self.cutOffT=t
		return self.cutOffT
	
	def setMaxfun(self,m):
		self.maxfun=m
		return self.maxfun
	
	def setOptTol(self,m):
		self.optTol=m
		return self.optTol
	
	def setLBD(self,b):
		self.LBD=b
		return self.LBD
	
	def setLBProd(self,b):
		self.LBProd=b
		return self.LBProd
	
	def setLBDegr(self,b):
		self.LBDegr=b
		return self.LBDegr
	
	def setUBD(self,b):
		self.UBD=b
		return self.UBD
	
	def setUBProd(self,b):
		self.UBProd=b
		return self.UBProd
	
	def setUBDegr(self,b):
		self.UBDegr=b
		return self.UBDegr
			
	def getEqu(self):
		return self.equOn
	
	def getFitPinned(self):
		return self.fitPinned
	
	def getFitProd(self):
		return self.fitProd
	
	def getFitDegr(self):
		return self.fitDegr
	
	def getSaveTrack(self):
		return self.saveTrack
	
	def getFitCutOffT(self):
		return self.fitCutOffT
	
	def getCutOffT(self):
		return self.cutOffT
	
	def getMaxfun(self):
		return self.maxfun
	
	def getOptTol(self):
		return self.optTol
	
	def getLBD(self):
		return self.LBD
	
	def getLBProd(self):
		return self.LBProd
	
	def getLBDegr(self):
		return self.LBDegr
	
	def getUBD(self):
		return self.UBD
	
	def getUBProd(self):
		return self.UBProd
	
	def getUBDegr(self):
		return self.UBDegr
	
	def setKineticTimeScale(self,s):
		self.kineticTimeScale=s
		return self.kineticTimeScale
	
	def getKineticTimeScale(self):
		return self.kineticTimeScale
	
	def setName(self,s):
		self.name=s
		return self.name
	
	def getName(self):
		return self.name
		
	def setX0D(self,x):
		self.x0[0]=x
		return self.x0[0]
	
	def setX0Prod(self,x):
		self.x0[1]=x
		return self.x0[1]
	
	def setX0Degr(self,x):
		self.x0[2]=x
		return self.x0[2]
	
	def getX0D(self):
		return self.x0[0]
	
	def getX0Prod(self):
		return self.x0[1]
	
	def getX0Degr(self):
		return self.x0[2]
	
	def setX0(self,x):
		self.x0=x
		return self.x0
	
	def checkPinned(self):
		b=True
		for i,r in enumerate(self.ROIsFitted):
			b = b + len(self.embryo.tvecData)==len(r.dataVecPinned) + len(self.embryo.simulation.tvecSim)==len(r.simVecPinned) 
		return b
	
	def checkSimulated(self):
		b=True
		for r in self.ROIsFitted:
			b = b + len(r.simVec)==len(self.embryo.simulation.tvecSim)
		return b
		
		
	def updateVersion(self):
		fittemp=fit(self.embryo,"temp")
		pyfrp_misc_module.updateObj(fittemp,self)
		return self