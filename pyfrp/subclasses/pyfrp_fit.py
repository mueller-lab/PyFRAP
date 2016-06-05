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
from pyfrp.modules import pyfrp_stats_module


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
		self.kineticTimeScale=1.
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
		self.RsqByROI={}
		
		#Equalization
		self.equFacts=[]
		
		#Empty result dataseries
		self.tvecFit=embryo.tvecData

	def addROI(self,r):
		
		"""Adds ROI to the list of fitted ROIs.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to be used for fitting.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		if r not in self.ROIsFitted:
			self.ROIsFitted.append(r)
		return self.ROIsFitted
	
	def addROIByName(self,name):
		
		"""Adds ROI to the list of fitted ROIs, given a specific name.
		
		Args:
			name (str): Name of ROI to be used for fitting.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		r=self.embryo.getROIByName(name)
		return self.addROI(r)
	
	def addROIById(self,Id):
		
		"""Adds ROI to the list of fitted ROIs, given a specific ROI Id.
		
		Args:
			Id (int): Id of ROI to be used for fitting.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		r=self.embryo.getROIById(Id)
		return self.addROI(r)
	
	def removeROI(self,r):
		
		"""Removes ROI from the list of fitted ROIs.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to be removed.
		
		Returns:
			list: Updated list of ROIs used for fitting.
			
		"""
		
		if r in self.ROIsFitted:
			self.ROIsFitted.remove(r)
		return self.ROIsFitted
	
	def getX0(self):
		
		"""Returns initial guess of fit in the form that is useful for 
		the call of the optimization algorithm.
		
		Copies x0 into local variable to pass to solver, pop entries that are 
		currently not needed since they are turned off via ``fitProd`` or ``fitDegr``.
		
		.. note:: Always gets executed at the start of ``run``.
		
		Returns:
			list: Currently used x0.
			
		"""
		
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
		
		"""Generates tuple of boundary tuples, limiting parameters 
		varied during SSD minimization.
		
		Will generate exactly the boundary tuple that is currently 
		useful to the optimization algorithm, meaning that only 
		values that are needed  since they are turned on via ``fitProd`` or ``fitDegr``
		will be included into tuple.
		
		Will use values that are stored in ``LBx`` and ``UBx``, where ``x`` is 
		``D``, ``Prod``, or ``Degr`` for the creation of the tuples. 
		
		.. note:: Always gets executed at the start of ``run``.
		
		Returns:
			tuple: Boundary value tuple.
			
		"""
	
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
		
	def run(self,debug=False,ax=None):
		
		"""Runs fit.
		
		Fitting is done by passing fit object to :py:func:`pyfrp.modules.pyfrp_fit_module.FRAPFitting`.
		This function then calls all necessary methods of fit to prepare it for optimization and
		then passes it to optimization algorithm.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			ax (matplotlib.axes): Axes to show debugging plots in.
		
		Returns:
			pyfrp.subclasses.pyfrp_fit.fit: ``self``.
		
		"""
		
		self=pyfrp_fit_module.FRAPFitting(self,debug=debug,ax=ax)
		return self
	
	def assignOptParms(self,res):
		
		r"""Assigns optimal parameters found by optimization algorithm to
		attributes in fit object depending on fit options chosen.
		
		Args:
			res (list): Result array from optimization algorithm.
		
		Returns:
			tuple: Tuple containing:
			
				* DOptPx (float): Optimal diffusion coefficient in :math:`\frac{\mathrm{px}^2}{s}}`.
				* prod (float): Optimal production rate in :math:`\frac{\[c\]}{s}}`.
				* degr (float): Optimal degradation rate in :math:`\frac{1}{s}}`.
				* DOptMu (float): Optimal diffusion coefficient in :math:`\frac{\mu\mathrm{m}^2}{s}}`.
				
		"""
		
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
	
	def plotFit(self,ax=None,legend=True,title=None):
		
		"""Plots fit, showing the result for all fitted ROIs.
		
		.. note:: If no ``ax`` is given, will create new one.
		
		.. image:: ../imgs/pyfrp_fit/fit.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes used for plotting.
			legend (bool): Show legend.
			title (str): Title of plot.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIsFitted:
			ax=r.plotFit(self,ax=ax,legend=legend,title=title)
			
		return ax
	
	def printResults(self):
		
		"""Prints out main results of fit."""
		
		printObjAttr('DOptMu',self)
		printObjAttr('DOptPx',self)
		printObjAttr('prodOpt',self)
		printObjAttr('degrOpt',self)
		printObjAttr('success',self)
		printObjAttr('Rsq',self)
		printObjAttr('MeanRsq',self)
		
		return True
	
	def setOptMeth(self,m):
		
		"""Sets optimization method.
		
		Available optimization methods are:
		
			* Constrained Nelder-Mead
			* Nelder-Mead
			* TNC
			* L-BFGS-B
			* SLSQP
			* brute
			* BFGS
			* CG
		
		See also http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.optimize.minimize.html and
		http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.optimize.brute.html#scipy.optimize.brute .
		
		You can find out more about the constrained Nelder-Mead algorithm in the documentation of 
		:py:func:`pyfrp.modules.pyfrp_optimization_module.constrObjFunc`.
		
		Args:
			m (str): New method.
			
		"""
		
		self.optMeth=m
		return self.optMeth
	
	def getOptMeth(self):
		
		"""Returns the currently used optimization algorithm.
		
		Returns:
			str: Optimization algorithm.
			
		"""
		
		return self.optMeth
		
	def isFitted(self):
		
		"""Checks if fit already has been run and succeeded.
		
		Returns:
			bool: ``True`` if success.
			
		"""
		
		return self.DOptMu!=None
		
	def setEqu(self,b):
		
		"""Turns on/off equalization.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.equOn=b
		return self.equOn
	
	def setFitPinned(self,b):
		
		"""Turns on/off if pinned series are supposed to be fitted.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.fitPinned=b
		return self.fitPinned
	
	def setFitProd(self,b):
		
		"""Turns on/off if production is supposed to be considered in fit.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.fitProd=b
		return self.fitProd
	
	def setFitDegr(self,b):
		
		"""Turns on/off if degradation is supposed to be considered in fit.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.fitDegr=b
		return self.fitDegr
	
	def setSaveTrack(self,b):
		
		"""Turns on/off if fitting process is supposed to be stored.
		
		This then can then be used to following the convergence of 
		the optimization algorithm and possibly to identify local minima.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		self.saveTrack=b
		return self.saveTrack
	
	def setFitCutOffT(self,b):
		
		"""Turns on/off if only a certain fraction of the timeseries
		is supposed to be fitted.
		
		.. warning:: This option is currently VERY experimental. Fitting might
		   crash.
		
		Args:
			b (bool): New flag value.
			
		Returns:
			bool: New flag value.
		
		"""
		
		
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
	
	def computeStats(self):
		self=pyfrp_stats_module.computeFitRsq(self)
		
	def printRsqByROI(self):
		print "Rsq Values by ROI for fit ", self.name
		printDict(self.RsqByROI)
		
		