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

"""Statistics module for PyFRAP toolbox, mainly used to evaluate goodness of fit, but also providing functions
to assess overall measurement statistics.
"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

import numpy as np

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def computeFitRsq(fit):
	
	"""Computes R-squared values for fit object.
	
	R-squared values contain:
		
		* Mean R-squared value over all ROIs included in fit, stored in ``fit.MeanRsq``.
		* Product of R-squared value over all ROIs included in fit, stored in ``fit.Rsq``.
		* R-squared value for each ROIs included in fit, stored in ``fit.RsqBuROI``.
		
	Args:
		fit (pyfrp.subclasses.pyfrp_fit.fit): Fit object.
		
	Returns:
		pyfrp.subclasses.pyfrp_fit.fit: Updated fit object.
	
	"""
	
	Rsqs=[]
	
	#Compute Rsqs for regions used for fitting	
	for i in range(len(fit.dataVecsFitted)):
		r=Rsq(fit.dataVecsFitted[i],fit.fittedVecs[i])
		fit.RsqByROI[fit.ROIsFitted[i].name]=r
		Rsqs.append(r)
	
	fit.MeanRsq=np.mean(Rsqs)	
	fit.Rsq=np.prod(Rsqs)
		
	return fit

def Rsq(data,x):
	
	r"""Computes R-squared values for fit series to data series.
	
	R-squared value is being computed as 
	
	.. math:: R^2 = 1 - \frac{\sum\limits_i (x_i - d_i)^2}{\sum\limits_i (d_i - \bar{d} )^2}
	
	Args:
		x (numpy.ndarray) Fit series.
		data (numpy.ndarray): Data series.
		
	Returns:
		float: R-squared value.
	
	"""
	
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

def computeSSD(data,x):
	
	r"""Computes sum of squared differences (SSD) of fit series to data series.
	
	The SSD is computed by
	
	.. math:: SSD = \sum\limits_i (x_i - d_i)^2
	
	Args:
		x (numpy.ndarray) Fit series.
		data (numpy.ndarray): Data series.
		
	Returns:
		float: SSD.
	"""
	
	return sum((data-x)**2)

def parameterStats(x):
	
	r"""Returns mean, standard deviation and standard error of array.
	
	Note that standard error is computed as 
	
	.. math:: \sigma_n = \frac{\sigma}{\sqrt{n}}
	
	where :math:`n` is the number of samples in ``x`` and :math: :math:`\sigma`
	is the standard deviation over ``x``.
	
	Args:
		x (numpy.ndarray): List of values.
	
	Returns:
		tuple: Tuple containing:
			
			* xMean (float): Mean of x.
			* xStd (float): Standard deviation of x.
			* xSterr (float): Standard error of x.
			
	"""
	
	return np.mean(x), np.std(x), np.std(x)/np.sqrt(len(x))
	


