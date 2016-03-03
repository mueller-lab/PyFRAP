#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Parameter fitting module for PyFRAP toolbox, including following functions:

#(1)  parm_fitting: Checks parameter fitting parameters, calls optimization algorithm and writes results into embryo object 
#(2)  calc_ssd: Objective function for minimization of SSD. Includes equalization steps, debugging plots and cut_off possibility. 
#(3)  pin_conc: Pins simulation results and dataseries between 0 and 1.

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import sys

#Numpy/Scipy
import numpy as np
from fipy import *
from scipy import interpolate
import scipy.optimize as sciopt

#PyFRAP
import pyfrp_stats_module
import pyfrp_plot_module 
import pyfrp_optimization_module 

from pyfrp_term_module import *


#matplotlib
import matplotlib.pyplot as plt

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fits scaling solution to data

def FRAPFitting(fit,debug=False,gui=None,ax=None):
		
	#Counter for function calls
	global iterations
	iterations=0

	#Building x0
	x0=fit.getX0()
	
	#Building bounds
	bnds=fit.getBounds()
		
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Calling optimization algorithms
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#Calling optimizers
	if fit.optMeth=='brute':
		res=optimize.brute(FRAPObjFunc, bnds,args=(fit,debug,ax,False), full_output=bool(debug),finish=optimize.fmin)
		
	elif fit.optMeth=='Constrained Nelder-Mead':
		
		x0=pyfrp_optimization_module.transformX0(fit.x0,[fit.LBD,fit.LBProd,fit.LBDegr],[fit.UBD,fit.UBProd,fit.UBDegr])
		res=sciopt.fmin(pyfrp_optimization_module.constrObjFunc,x0,args=(fit,debug,ax,False),ftol=fit.optTol,maxiter=fit.maxfun,disp=bool(debug),full_output=True)
	
	elif fit.optMeth=='Anneal':
		random.seed(555)
		res=sciopt.minimize(FRAPObjFunc, x0,args=(fit,debug,ax,False), method='Anneal')
	else:
		res=sciopt.minimize(FRAPObjFunc,x0,args=(fit,debug,ax,False),method=fit.optMeth,tol=fit.optTol,options={'maxiter': fit.maxfun, 'disp': bool(debug)})
	
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Run for one last time to get final fit
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	if fit.optMeth=='Constrained Nelder-Mead':
		resNew=pyfrp_optimization_module.xTransform(res[0],[fit.LBD,fit.LBProd,fit.LBDegr],[fit.UBD,fit.UBProd,fit.UBDegr])
		
		fit=FRAPObjFunc(resNew,fit,debug,ax,True)
		
	elif fit.optMeth=='brute':
		fit=FRAPObjFunc(res[0],fit,debug,ax,True)
		
	else:	
		fit=FRAPObjFunc(res.x,fit,debug,ax,True)
	
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Saving results in embryo object
	#------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
	if fit.optMeth=='brute':
		fit.assignOptParms(res[0])
			
		fit.SSD=res[1]
		fit.success=True
		
		fit.iterations=iterations
		#In bruteforce iterations = fcalls???
		fit.fcalls=iterations
		
	elif fit.optMeth=='Constrained Nelder-Mead':
		
		resNew=pyfrp_optimization_module.xTransform(res[0],[fit.LBD,fit.LBProd,fit.LBDegr],[fit.UBD,fit.UBProd,fit.UBDegr])
		
		fit.assignOptParms(resNew)
		
		fit.SSD=res[1]
		fit.success=not bool(res[4])
		fit.iterations=res[2]
		fit.fcalls=res[3]
		
	else:	
		fit.assignOptParms(res.x)
			
		fit.SSD=res.fun
		fit.success=res.success
		
		fit.iterations=res.nit
		fit.fcalls=res.nfev
		
	fit=pyfrp_stats_module.computeFitRsq(fit)

	return fit

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def scaleTime(tvec,D,Dnew):
	return D/Dnew*tvec

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def interpolateSolution(tvecData,tvecScaled,yvec):
	fscal=interpolate.interp1d(tvecScaled,yvec)
	yvecScaled=fscal(tvecData)
	return yvecScaled

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def assignInputVariables(x,fit):
	Dnew=x[0]
	if fit.fitProd and fit.fitDegr:
		prod=x[1]
		degr=x[2]
	
	elif fit.fitProd and not fit.fitDegr:
		prod=x[1]
		degr=fit.x0[2]
			
	elif not fit.fitProd and fit.fitDegr:
		degr=x[1]
		prod=fit.x0[1]
		
	elif not fit.fitProd and not fit.fitDegr:
		degr=fit.x0[2]
		prod=fit.x0[1]
	
	return Dnew,prod,degr

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def checkInput(x,iteration,fit):
	x=np.array([1,fit.fitProd,fit.fitDegr])*x
	if min(x)<0:
		if iteration==0:
			#If initial guess is out of bounds, just do nothing
			print printWarning("Check your initial guess, some value is negative.")
			return False
		else:
			print printWarning("One of the values of x is negative")
			return False
	return True
		
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def downscaleKinetics(prod,degr,rate):
	return float(prod)/rate, float(degr)/rate

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def getTvecCutIndex(tvec,tCut):
	return max(np.where(tvec<tCut)[0])

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def cutTvec(tvec,tCutIndex):
	return tvec[:tCutIndex]

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def scaleROIs(fit,Dnew):
	
	scaledSimVecs=[]
	dataVecs=[]
		
	#Rescaling the time series by D/Dnew
	tvecScaled=scaleTime(fit.embryo.simulation.tvecSim,fit.embryo.simulation.D,Dnew)
	tvecData=fit.embryo.tvecData
	
	#If cutOff option is given, need to find time point and corresponding steps of cutOff
	if fit.fitCutOffT:
		fit.cutOffStepData=getTvecCutIndex(fit.embryo.tvecData)
		fit.cutOffStepSim=getTvecCutIndex(fit.embryo.simulation.tvecSim)
		tvecScaled=scaleTime(fit.embryo.simulation.tvecSim[:fit.cutOffStepSim],fit.embryo.simulation.D,Dnew)
		tvecData=np.array(fit.embryo.tvecData[:fit.cuttOffStepData])
			
	#Interpolating new solution		
	for r in fit.ROIsFitted:
		if fit.fitPinned:
			if fit.fitCutOffT:
				scaledSimVecs.append(np.array(interpolateSolution(fit.embryo.tvecData,tvecScaled,r.simVecPinned[:fit.cutOffStepSim])))
				dataVecs.append(np.array(r.dataVecPinned[:fit.cutOffStepData]))
			else:
				scaledSimVecs.append(np.array(interpolateSolution(fit.embryo.tvecData,tvecScaled,r.simVecPinned)))
				dataVecs.append(np.array(r.dataVecPinned[:fit.cutOffStepData]))
		else:
			if fit.fitCutOffT:
				scaledSimVecs.append(interpolateSolution(np.array(fit.embryo.tvecData,tvecScaled,r.simVec[:fit.cutOffStepSim])))
				dataVecs.append(np.array(r.dataVecPinned))
			else:
				scaledSimVecs.append(np.array(interpolateSolution(fit.embryo.tvecData,tvecScaled,r.simVec)))
				dataVecs.append(np.array(r.dataVec))
	
	return fit,tvecScaled,tvecData,scaledSimVecs,dataVecs

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def addKineticsToSolution(scaledSimVecs,tvec,prod,degr):
	
	rescaledSimVecs=[]
	
	for scaledSimVec in scaledSimVecs:
	
		#Both production and degredation
		if prod>0 and degr>0:
			rescaledSimVec=scaledSimVec*np.exp(-degr*tvec)-(prod/degr)*np.exp(-degr*tvec)+(prod/degr)
		
		#Just production
		elif prod>0 and degr<=0:
			rescaledSimVec=scaledSimVec+prod*tvec
			
		#Just degradation
		elif prod<=0 and degr>0:
			rescaledSimVec=scaledSimVec*np.exp(-degr*tvec)
		
		#Neither production nor degradation
		rescaledSimVec=scaledSimVec
		
		rescaledSimVecs.append(rescaledSimVec)
		
	return rescaledSimVecs	

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def computeEquFactors(dataVec,simVec):
			
	#Add little Noise to concentration profiles so we don't get singularities
	simVec=simVec+1E-10
	
	#Calculate equalization factors
	equFacts=simVec/(dataVec+1E-10)
	
	""" 
	Check requirements:
	these two if statements set the lower and upper bounds for the equalization factor.
	In practical terms, this means that the volume fraction is not allowed to increase by
	more than 3fold (>3), and that the immobile fraction cannot be more than 90% (< 0.1)
	"""
	
	equFacts[np.where((equFacts>3.) + (equFacts<0.1))[0]]=1.
	
	return equFacts

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def equalize(dataVecs,simVecs):
	
	equSimVecs=[]
	equFacts=[]
	
	#Compute Equalization factors
	for i in range(len(dataVecs)):
		equFacts.append(computeEquFactors(dataVecs[i],simVecs[i]))
	
	#Compute SSD for all ROIs used for fitting for all equalization factors
	SSDs=[]
	for j in range(len(equFacts[0])):
		ssds=[]
		for i in range(len(equFacts)):	
			equSimVec=simVecs[i]/equFacts[i][j]
			ssds.append(pyfrp_stats_module.computeSSD(dataVecs[i],equSimVec))
			
		SSDs.append(sum(ssds))
	
	#Compute Final Equalization Factor (the one that minimizes SSD)
	equFactFinalIdx=SSDs.index(min(SSDs))
	
	#Return equalized simulation vectors
	for i,vec in enumerate(simVecs):
		equSimVecs.append(vec/equFacts[i][equFactFinalIdx])
		
	return equSimVecs
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Objective function for fitting
	
def FRAPObjFunc(x,fit,debug,ax,returnFit):
	
	#Counting function calls
	global iterations
	iterations=iterations+1

	#Check if any variable is negative
	if not checkInput(x,iterations,fit):
		return 2*fit.SSD
	
	#Assign Input Values
	Dnew,prod,degr = assignInputVariables(x,fit)

	#Rescaling degr and prod
	prod,degr = downscaleKinetics(prod,degr,fit.kineticTimeScale)
	
	if debug:
		print "------------------------------------------"
		print "Dnew=",Dnew, "prod=", int(fit.fitProd)*prod, "degr=", int(fit.fitDegr)*degr
	
	#Scale simulation vectors
	try:
		fit,tvecScaled,tvecData,scaledSimVecs,dataVecs = scaleROIs(fit,Dnew)
	except ValueError:
		printWarning("Scaling failed with Dnew = " + str(Dnew))
		if returnFit:
			return fit
		else:	
			return 100000000
	
	#Add Kinetics
	scaledSimVecs =  addKineticsToSolution(scaledSimVecs,tvecData,prod,degr)
	
	#Equalize
	if fit.equOn:
		scaledSimVecs=equalize(dataVecs,scaledSimVecs)
	
	#Compute final SSD
	ssds=[]
	for i in range(len(dataVecs)):
		ssds.append(pyfrp_stats_module.computeSSD(dataVecs[i],scaledSimVecs[i]))
	SSD=sum(ssds)
	
	#Write final fitting vectors in ROI objects
	fit.fittedVecs=scaledSimVecs
	fit.dataVecsFitted=dataVecs
	
	#Live-Plot
	if debug:
		if iterations==1:
			if ax==None:
				global fig
				global axMon
				fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Live Fitting Output"],sup=fit.name)
				axMon=axes[0]
			else:
				axMon=ax
				
		else:
			axMon.cla()
		
		for i,r in enumerate(fit.ROIsFitted):
			axMon.plot(tvecData,dataVecs[i],color=r.color,linestyle='-')
			axMon.plot(tvecData,scaledSimVecs[i],color=r.color,linestyle='--')
			
		plt.draw()
		plt.pause(0.00001)
	
	fit.SSD=SSD
	fit.tvecFit=tvecData
	
	if returnFit:
		return fit
	else:
		return SSD

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def computePinVals(vec,useMin=False,useMax=False,bkgdVal=None,debug=False):
	if bkgdVal==None:
		bkgdVal = computeBkgd(vec,useMin=useMin,debug=debug)
	vec2=np.asarray(vec)-bkgdVal
	normVal = computeNorm(vec2,useMax=useMax,debug=debug)
	
	return bkgdVal,normVal

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
def computeBkgd(vec,useMin=False,debug=False):
	minVal=min(vec)
	if debug and not useMin:
		if minVal<vec[0]:
			printWarning("First value of vector is not the smallest. You may want to use useMin=True to get an optimal background value.")
	if useMin:
		return minVal
	else:
		return vec[0]

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def computeNorm(vec,useMax=False,debug=False):
	maxVal=max(vec)
	if debug and not useMax:
		if maxVal<vec[-1]:
			printWarning("Last value of vector is not the maximum. You may want to use useMax=True to get an optimal background value.")
	if useMax:
		return maxVal
	else:
		return vec[-1]
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Substract background and normalize: Pin concentrations between 0 and 1:

def pinConc(vec,bkgdVal,normVal,switch_thresh=0.95,axes=None,debug=False,tvec=None,color='b'):
	
	#Copying unpinned timeseries
	vecPinned=list(vec)
		
	#Substracting background from all timeseries
	vecPinned[:]=[x-bkgdVal for x in vecPinned]
	
	if debug:
		
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,2],titles=['After substraction','After norming'],sup="pinConc Debug Plot")
			
		if tvec==None:
			tvec=np.arange(np.shape(vecPinned))
		
		axes[0].plot(tvec,vecPinned,color=color,linestyle='-',label="pinned vector") 
		axes[0].plot(tvec,vec,color=color,linestyle='--',label="data vector") 
		
		#axes[0].set_ylim([0,1.5])
		
		plt.draw()
		
	#Dividing background from all timeseries
	vecPinnedOld=list(vecPinned)
	vecPinned[:]=[x/normVal for x in vecPinned]
	
	if debug:
		
		axes[1].plot(tvec,vecPinnedOld,color=color,linestyle='-',label="pinned vector") 
		axes[1].plot(tvec,vecPinned,color=color,linestyle='-.',label="normed vector") 
		
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		#pyfrp_plot_module.redraw(axes[0])
		plt.draw()
		
	return vecPinned





	