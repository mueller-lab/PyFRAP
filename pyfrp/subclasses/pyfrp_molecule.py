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

#Main molecule class saving all important data, parameters and results for particular molecule.
	
#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_IO_module
from pyfrp.modules import pyfrp_stats_module

#PyFRAP Classes
import pyfrp_embryo

#Standard packages
import os


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Molecule object

class molecule:
	
	#Creates new molecule object
	def __init__(self,name):
		
		#Data
		self.name=name
		self.embryos=[]
		
		#Results
		self.selFits=[]
		self.DOptMu=None
		self.DOptPx=None
		self.DOptPxStd=None
		self.DOptPxSterr=None
		self.prodOpt=None
		self.degrOpt=None
		self.MeanRsq=None
		self.Rsq=None
		self.DMuStd=None
		self.prodStd=None
		self.degrStd=None
		self.DMuStErr=None
		self.prodStErr=None
		self.degrStErr=None
		
		self.crucialParameters=["equOn","fitPinned","fitProd","fitDegr","LBD","LBProd","LBDegr","UBD","UBProd","UBDegr"]
			
		#self.fitting_parms=[]
	
		
				
	def addEmbryo(self,embryo):
		self.embryos.append(embryo)
		return self.embryos
	
	def newEmbryo(self,name):
		emb=pyfrp_embryo.embryo(name)
		self.embryos.append(emb)
		return self.embryos[-1]
		
	def removeEmbryo(self,embryo):
		self.embryos.remove(embryo)
		return self.embryos		
	
	def save(self,fn=None):
		
		if fn==None:
			fn=self.name+".pk"
		
		pyfrp_IO_module.saveToPickle(self,fn=fn)	
		
		return fn
	
	def extractEmbryos2Files(self,fn=""):
		
		if fn=="":
			try:
				os.system("mkdir -v embryoFiles/")
				fn="embryoFiles/"
			except OSError:
				pass
		
		for embryo in self.embryos:
			embryo.save(fn+embryo.getName()+".emb")
		
	def updateVersion(self):
		
		#Create temporarly a blank molecule file
		moltemp=molecule("temp")
		
		#Update molecule file
		pyfrp_misc_module.updateObj(moltemp,self)
		
		#Update all embryos and fits
		for emb in self.embryos:
			
			emb.updateVersion()
			
			for fit in emb.fits:
				fit.updateVersion()
		
		return self
		
	def getEmbryoByName(self,s):
		
		for emb in self.embryos:
			if str(s) == emb.name:
				return emb
		return False
	
	def setName(self,n):
		self.name=n
		return self.name
	
	def getName(self):
		return self.name
	
	def sumUpResults(self,sameSettings=False):
		
		if sameSettings:
			
			lastFit=self.selFits[0]
			
			#Check that all fits have roughly the same parameters
			for fit in self.selFits:
				
				#Check if fit is fitted
				if not fit.isFitted():
					printError("Cannot average fits, since fits " + fit.name + "has not been fitted yet.")
					return False
				
				same,different,notInBoth=compareObjAttr(lastFit,fit)
				
				for item in self.crucialParameters:   
					if item in different.keys():
						printError("Cannot average fits, since fits " + lastFit.name + " and " + fit.name + "do not have the same value for " + item +". However, this parameter is marked as crucial.")
						return False
					elif item in notInBoth:
						printError("Cannot average fits, since one of the fits " + lastFit.name + " and " + fit.name + " lacks the parameter " + item +".")
						return False
		
		#Compute Statistics
		self.DOptMu,self.DOptMuStd,self.DOptMuSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"DOptMu"))
		self.DOptPx,self.DOptPxStd,self.DOptPxSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"DOptPx"))
		self.prodOpt,self.prodOptStd,self.prodOptSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"prodOpt"))
		self.degrOpt,self.degrOptStd,self.degrOptSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"degrOpt"))
		
		self.Rsq=np.mean(pyfrp_misc_module.objAttrToList(self.selFits,"Rsq"))
		self.MeanRsq=np.mean(pyfrp_misc_module.objAttrToList(self.selFits,"MeanRsq"))
			
		return	True
		
		