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

#Statistics module for PyFRAP toolbox, including following functions:
#(1) fit_Rsq: Calulates Rsq value of fit

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================


import numpy as np

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calulates Rsq value of fit
	
def computeFitRsq(fit):
	Rsqs=[]
	
	#Compute Rsqs for regions used for fitting	
	for i in range(len(fit.dataVecsFitted)):
		r=Rsq(fit.dataVecsFitted[i],fit.fittedVecs[i])
		fit.RsqByROI[fit.ROIsFitted[i].name]=r
		Rsqs.append(r)
	
	fit.MeanRsq=np.mean(Rsqs)	
	fit.Rsq=np.prod(Rsqs)
		
	return fit

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calulates Rsq value of fit

def Rsq(data,x):
	
	#Make numpy array to avoid list substraction problem
	data=np.asarray(data)
	x=np.asarray(x)
	
	#Compute ssd between x and data
	ssd=computeSSD(data,x)
	
	#Compute mean data
	meanData=np.mean(data)
	
	#Get derivation of datapoints from mean
	SStot=sum((data-meanData)**2)
	
	#Calculate RSq value
	Rsq=1-ssd/SStot
	
	return Rsq

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calulates Rsq value of fit

def computeSSD(data,x):
	return sum((data-x)**2)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns mean/std/stderr of array

def parameterStats(x):
	return np.mean(x), np.std(x), np.std(x)/np.sqrt(len(x))
	


