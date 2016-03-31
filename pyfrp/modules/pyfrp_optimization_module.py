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

#Parameter fitting module for PyFRAP toolbox, including following functions:

#(1)  parm_fitting: Checks parameter fitting parameters, calls optimization algorithm and writes results into embryo object 


#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP
import pyfrp_fit_module 

from pyfrp_term_module import *

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Substract background and normalize: Pin concentrations between 0 and 1:

def constrObjFunc(x,fit,debug,ax,returnFit):
	
	x=xTransform(x,[fit.LBD,fit.LBProd,fit.LBDegr],[fit.UBD,fit.UBProd,fit.UBDegr])

	ssd=pyfrp_fit_module.FRAPObjFunc(x,fit,debug,ax,returnFit)
	
	return ssd

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Transforms unconstrained problem into constrained problem
	
def xTransform(x,LB,UB):
	
	#Make sure everything is float
	x=np.asarray(x,dtype=np.float64)
	LB=np.asarray(LB,dtype=np.float64)
	UB=np.asarray(UB,dtype=np.float64)
	
	#Check if LB_D==0, then add a little noise to it so we do not end up with xtrans[D]==0 and later have singularities when scaling tvec
	if LB[0]==0:
		LB[0]=1E-10
	
	#Determine number of parameters to be fitted
	nparams=len(x)

	#Make empty vector
	xtrans = np.zeros(np.shape(x))
	
	# k allows some variables to be fixed, thus dropped from the
	# optimization.
	k=0

	for i in range(nparams):

		#Upper bound only
		if UB[i]!=None and LB[i]==None:
		
			xtrans[i]=UB[i]-x[k]**2
			k=k+1
			
		#Lower bound only	
		elif UB[i]==None and LB[i]!=None:
			
			xtrans[i]=LB[i]+x[k]**2
			k=k+1
		
		#Both bounds
		elif UB[i]!=None and LB[i]!=None:
			
			xtrans[i] = (np.sin(x[k])+1.)/2.*(UB[i] - LB[i]) + LB[i]
			xtrans[i] = max([LB[i],min([UB[i],xtrans[i]])])
			k=k+1
		
		#No bounds
		elif UB[i]==None and LB[i]==None:
		
			xtrans[i] = x[k]
			k=k+1
			
		#Note: The original file has here another case for fixed variable, but since we made the decision earlier which when we call frap_fitting, we don't need this here.
	
	return xtrans	

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#transform x0 from unconstrained to constrained problem
		
def transformX0(x0,LB,UB):

	# transform starting values into their unconstrained surrogates. Check for infeasible starting guesses.
	x0u = list(x0)
	
	nparams=len(x0)
	
	k=0
	for i in range(nparams):
		
		#Upper bound only
		if UB[i]!=None and LB[i]==None:
			if UB[i]<=x0[i]:
				x0u[k]=0
			else:
				x0u[k]=sqrt(UB[i]-x0[i])	
			k=k+1
			
		#Lower bound only
		elif UB[i]==None and LB[i]!=None:
			if LB[i]>=x0[i]:
				x0u[k]=0
			else:
				x0u[k]=np.sqrt(x0[i]-LB[i])	
			k=k+1
		
		
		#Both bounds
		elif UB[i]!=None and LB[i]!=None:
			if UB[i]<=x0[i]:
				x0u[k]=np.pi/2
			elif LB[i]>=x0[i]:
				x0u[k]=-np.pi/2
			else:
				x0u[k] = 2*(x0[i] - LB[i])/(UB[i]-LB[i]) - 1;
				#shift by 2*pi to avoid problems at zero in fminsearch otherwise, the initial simplex is vanishingly small
				x0u[k] = 2*np.pi+np.arcsin(max([-1,min(1,x0u[k])]));
			k=k+1
		
		#No bounds
		elif UB[i]==None and LB[i]==None:
			x0u[k] = x[i]
			k=k+1
	
	return x0u