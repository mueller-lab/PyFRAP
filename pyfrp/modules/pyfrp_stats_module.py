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
	


