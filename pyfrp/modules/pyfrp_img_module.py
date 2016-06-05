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

"""Image analysis module for PyFRAP toolbox. This is one of the key modules of PyFRAP, gathering all necessary 
image manipulation and reading functions used by PyFRAP. Mainly focuses on:

	* Analyzing complete FRAP datasets as specified in a :py:class:`pyfrp.subclasses.pyfrp_analysis.analysis` instance.
	* Reading concentrations from images over specific regions as specified in  :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` instances.
	* Handling and saving rim concentrations.
	* Image manipulation through filters/quadrant reduction.
	* Image plotting and histograms.
	* Fiji wrapping.
	
"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#numpy
import numpy as np

#Plotting
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Misc
import sys
import time
import os

#PyFRAP modules
import pyfrp_misc_module
import pyfrp_plot_module
import pyfrp_idx_module
from pyfrp_term_module import *

#Image processing
import skimage
import skimage.morphology
import skimage.io
#skimageVersion=skimage.__version__
try:
	if int(skimage.__version__.split('.')[1])<11:
		import skimage.filter as skifilt
	else:
		import skimage.filters  as skifilt
except:
	import skimage.filters as skifilt
try:	
	import scipy.signal as spsig
except:
	printWarning("Cannot import scipy.signal. Will not be able to do spsig.gaussian")

                    
#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def analyzeDataset(analysis,signal=None,embCount=None,debug=False,debugAll=False,showProgress=True):
	
	"""Main dataset analysis function doing the following steps.
	
	* Reset all data arrays for all ROIs.
	* Computes flattening/norm/background mask if necessary.
	* Loops through images, processing images and computing mean concentraions per ROI.
	* Showing final debugging plots if selected.
	
	Args:
		analysis (pyfrp.subclasses.pyfrp_analysis):  Object containing all necessary information for analysis.
		
	Keyword Args:
		signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
		embCount (int): Counter of counter process if multiple datasets are analyzed. 
		debug (bool): Print final debugging messages and show debugging plots.
		debugAll (bool): Print debugging messages and show debugging plots of each step.
		showProgress (bool): Print out progress.
		
	Returns:
		pyfrp.subclasses.pyfrp_analysis: Performed analysis.
	"""
	
	if debug or signal==None:
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "Starting analyzing dataset " + analysis.embryo.name
		if debug:
			print 'Anaylsis options:'
			printDict(analysis.process)
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	#Set back all timeseries vectors
	for r in analysis.embryo.ROIs:
		r.resetDataVec()
	
	#Compute flattening mask if needed
	if 'flatten' in analysis.process.keys():
		flatteningMask=analysis.computeFlatteningMask()
	else:
		flatteningMask=None
	
	#Compute background mask if needed
	if 'bkgd' in analysis.process.keys():
		bkgdMask = analysis.computeBkgdMask(flatteningMask)
	else:
		bkgdMask = None
	
	#Load preimage if needed
	if 'norm' in analysis.process.keys():
		preMask = analysis.computePreMask(flatteningMask,bkgdMask)
	else:
		preMask = None
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Loop through images and compute concentrations
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	for i in range(len(analysis.embryo.fileList)):
		
		#Compose filename
		fnImg=str(analysis.embryo.getDataFolder()+'/'+analysis.embryo.getFileList()[i])
			
		#Reading in image
		img=loadImg(fnImg,analysis.embryo.dataEnc)
		
		#Check if skimage reads in image as 2D array, if not grab channel of image with maximum range
		if len(np.shape(img))>2:
			img,ind_max=getMaxRangeChannel(img,debug=debugAll)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Process image / Get image ready for readout
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		#Process Image
		img = processImg(img,analysis.process,flatteningMask,bkgdMask,preMask,analysis.dataOffset,debug=debugAll)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Compute concentrations
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		#Get rim concentrations
		concRim=getRimConc(analysis.embryo.ROIs,img,debug=debugAll)
			
		#Get concentrations in all ROIs
		for r in analysis.embryo.ROIs:
			r.dataVec.append(meanExtConc(r.imgIdxX,r.imgIdxY,img,concRim,r.numExt,analysis.addRimImg,debug=debugAll))
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Save first image and its concRim for simulation
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		if i==0:
			
			if analysis.embryo.simulation!=None:
				if 'quad' in analysis.process.keys():
					analysis.embryo.simulation.ICimg=np.flipud(convSkio2NP(flipQuad(img)))
				else:
					analysis.embryo.simulation.ICimg=convSkio2NP(img)
				
			analysis.concRim=concRim	
			
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Print Progress
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		if showProgress:
			currPerc=int(100*i/float(len(analysis.embryo.fileList)))
			
			if signal==None:
				sys.stdout.write("\r%d%%" %currPerc)  
				sys.stdout.flush()
			else:	
				if embCount==None:
					signal.emit(currPerc)
				else:
					signal.emit(currPerc,embCount)
	print
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Final debugging plots
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plot_module.makeSubplot([1,2],titles=["Analyzed timeseries", "ICs"],sup="analayzeDataset debugging plot",tight=True)
			
		#Plot concentration timeseries
		analysis.embryo.plotAllData(ax=axes[1],legend=True)
		
		#Plot ICs		
		analysis.embryo.simulation.showICimg(ax=axes[0])
		analysis.embryo.showAllROIBoundaries(ax=axes[0])
		
		#Draw
		plt.draw()
	
	return analysis

def convSkio2NP(img):
	
	"""Returns mean concentration over given indices. 

	Args:
		idxX (list): x-indices of pixels used.
		idxY (list): y-indices of pixels used. 
		vals (numpy.ndarray): Input image.
		
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		float: Mean concentration.
	"""
	
	return np.asarray(img)
	
def genFakeIC(res,valIn,valOut,offset,sidelength,radius,center,rim,add_rim_from_radius,debug=False):
	
	"""Create fake image consisting of circular ROI with pixel value valOut
	and a square ROI centered inside circle with pixel value valIn. 
	"""
	
	#Empty picture
	vals=np.nan*np.ones((res,res))
	
	#getting indices for regions
	idxCircleX,idxCircleY = pyfrp_idx_module.getCircleIdxImg(center,radius,res,debug=debug)
	idxSquareX,idxSquareY = pyfrp_idx_module.getSquareIdxImg(offset,sidelength,res,debug=debug)
	
	#Assign values
	vals=pyfrp_idx_module.ind2mask(vals,idxSquareX,idxSquareY,valIn)
	vals=pyfrp_idx_module.ind2mask(vals,idxCircleX,idxCircleY,valOut)
			
	#Debugging plot	
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Fake IC"],sup="genFakeIC debugging output")
		
		#Contour plot of fake ICs
		contplot=axes[0].imshow(vals,cmap=plt.cm.jet)
		cbar = plt.colorbar(contplot)
		
		#Draw
		plt.draw()
	
	return vals


def meanConc(idxX,idxY,vals,debug=False):
	
	"""Returns mean concentration over given indices. 

	Args:
		idxX (list): x-indices of pixels used.
		idxY (list): y-indices of pixels used. 
		vals (numpy.ndarray): Input image.
		
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		float: Mean concentration.
	"""
	
	if debug:
		print "======= mean_conc debugging output ======="
		print "len(idxX)", len(idxX)
		print "vals.min(), vals.max(), mean(vals)",vals.min(), vals.max(), mean(vals)
		print "vals[idxX,idxY].min(), vals[idxX,idxY].max(), mean(vals[idxX,idxY]) ",vals[idxX,idxY].min(), vals[idxX,idxY].max(), mean(vals[idxX,idxY])
		
	return np.mean(vals[idxX,idxY])

def meanExtConc(idxX,idxY,img,concRim,numExt,addRimImg,debug=False):
	
	"""Returns mean concentration, taking potential pixels outside of image into account. 

	Args:
		idxX (list): x-indices of pixels used.
		idxY (list): y-indices of pixels used. 
		img (numpy.ndarray): Input image.
		concRim (float): Rim concentration applied to pixels outside of image.
		numExt (int): Number of pixels outside of image.
		addRimImg (bool): Add rim concentraion.
		
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		float: Mean extendended concentration.
	"""
	
	#Adding up concentrations of all selected indices (can't use mean conc since addRimImg could be True)
	concSum=img[idxX,idxY].sum()
	concNum=len(idxX)
	
	#Check if I want to add rim
	if addRimImg:	
		
		#We assume that the concentration in all extended pixels is equals the concRim and add it numExt times to overall concentration 
		concSum=concSum+numExt*concRim
		concNum=concNum+numExt
	
	#Return final concentration
	return float(concSum)/float(concNum)


def getRimConc(ROIs,img,debug=False):
	
	"""Computes mean rim concentration from ROIs that have *useForRim* flag on. 

	Args:
		ROIs (list): List of pyfrp.subclasses.pyfrp_ROI objects.
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		debug (bool): Print debugging messages and show debugging plots.
		
	Returns:
		float: Rim concentration.
	"""
	
	#Loop through all ROIs and find the ones needed to compute concentration rim
	idxX=[]
	idxY=[]
	for r in ROIs:
		if r.useForRim:
			
			#Fuse all indices
			idxX=idxX+list(r.imgIdxX)
			idxY=idxY+list(r.imgIdxY)
			
			#Remove doubles
			idxX,idxY=pyfrp_idx_module.remRepeatedImgIdxs(idxX,idxY,debug=debug)
			
	#Compute concentration 
	return meanConc(idxX,idxY,img,debug=debug)

def flipQuad(img,debug=False,testimg=False):
	
	"""Flip image into quaddrant.

	Args:
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		debug (bool): Print debugging messages and show debugging plots.
		testimg (bool): Check with test image if flip works properly
		
	Returns:
		numpy.ndarray: Flipped image.
	"""
	
	if debug:
		
		if testimg:
			
			#Generate fake ICs with optimal settings
			img=genFakeIC(img.shape[0],0,1,[256.5-50,256.5-50],100,256,[256.5,256.5],101,0)
			
			#fixed_thresh(img,)
			
			#Check symmetry
			if not symmetryTest(img,debug=debug):
				raw_input()
				
		#Grab max/min values of original image	
		org_max=img.max()
		org_min=img.min()
		
		#Save original image
		org_img=img.copy()
		
		#Create figure
		fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["Original image","Flipped image before normalization","Flipped image after normalization"],sup="flipQuad debugging output")
		a=axes[0].imshow(img,vmin=org_min,vmax=org_max)
		plt.colorbar(a)
			
	#Grab image resolution
	res=np.shape(img)[0]
		
	#Add left side onto rigth side
	img=np.fliplr(img[:,:res/2]) + img[:,res/2:]

	#Add bottom side onto upper side
	img=img[:res/2,:] + np.flipud(img[res/2:,:])
	
	if debug:
		a=axes[1].imshow(img)
		plt.colorbar(a)
		
	#Normalize
	img=img/4.
	
	if debug:
	
		a=axes[2].imshow(img)
		
		plt.draw()
		plt.colorbar(a)
		
		
		print "======= flipQuad debugging output ======="
		print "Corner check: "
		print "Original Image corners: ",  org_img[0,0],org_img[0,-1],org_img[-1,0],org_img[-1,-1]
		print "Original Image average: ", (org_img[0,0]+org_img[0,-1]+org_img[-1,0]+org_img[-1,-1])/4.
		print "Flipped Image corners: ", img[0,0],img[0,-1],img[-1,0],img[-1,-1]
		print "Flipped image: ", img[-1,-1]
		
		print
		
		print "Unique check: "
		print "Original Image #values ", unique(org_img)
		print "Flipped Image #values ", unique(img)
		
		
		raw_input()
		
	return img

def unflipQuad(img,debug=False,testimg=False):
	
	"""Unflip image from quaddrant into normal picture.

	Args:
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		debug (bool): Print debugging messages and show debugging plots.
		testimg (bool): Check with test image if unflip works properly
		
	Returns:
		numpy.ndarray: Unflipped image.
	"""
	
	if debug:
		if testimg:
			
			#Generate fake ICs with optimal settings
			img=gen_fake_IC(2*img.shape[0],0,1,[256.5-50,256.5-50],100,256,[256.5,256.5],101,0)
			
			#Flip fake ICs
			img=flipQuad(img)
				
		#Grab max/min values of original image	
		org_max=img.max()
		org_min=img.min()
		
		#Save original image
		org_img=img.copy()
		
		#Create figure
		fig,axes = pyfrp_plot_module.makeSubplot([2,2],titles=["Original image","Before first flip","After first flip","After second flip"],sup="unflipQuad debugging output")
	
		a=axes[0].imshow(img,vmin=org_min,vmax=org_max)
		plt.colorbar(a)
	
	#Make new empty image
	new_img=np.zeros((2*img.shape[0],2*img.shape[1]))
	
	#Grab image resolution
	res=np.shape(new_img)[0]
	
	#Assign first quadrant to empty image
	new_img[:res/2,res/2:]=img
	img=new_img
	
	#Plot result
	if debug:
		a=axes[1].imshow(img)
		
	#Flip image UD and add to itself
	img=img+np.flipud(img)
	
	#Plot result
	if debug:
		a=axes[2].imshow(img)
		
	#Flip image LR and add to itself
	img=img+np.fliplr(img)
	
	#Plot result
	if debug:
		a=axes[3].imshow(img)
		
		plt.draw()
		
	#Run symmetry test in the end
	if debug:
		symmetryTest(img,debug=debug)
	
	return img
	
def symmetryTest(img,debug=False):
	
	"""Checks if images is LR and UD symmetirc. If so, returns True.

	Args:
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		debug (bool): Print debugging messages and show debugging plots.
	
	Returns:
		bool: True if image is both LR and UD symmetric
	"""
	
	#Grab image resolution
	res=np.shape(img)[0]
	
	#Create figure
	if debug:
		fig,axes = pyfrp_plot_module.makeSubplot([1,2],titles=["LR check", "UD check"],sup="symmetry_test debugging output")
		
	#Check LR flip
	lr=np.fliplr(img[:,:res/2]) - img[:,res/2:]
			
	#Draw LR flip
	if debug:
		a=axes[0].imshow(lr)
		plt.colorbar(a)
		
	#Check UD flip
	ud=img[:res/2,:] - np.flipud(img[res/2:,:])

	#Draw UD flip
	if debug:
		a=axes[1].imshow(ud)
		plt.colorbar(a)
	
	#Draw
	plt.draw()
	
	if debug:
		if (lr.sum()+ud.sum())==0:
			pass
		else:
			print "======= symmetry_test debugging output ======="
			print "Image is not symmetric"
	
	#Return 1 if symmetric for both UD and LR
	
	return not bool(lr.sum()+ud.sum())

def getMaxRangeChannel(img,debug=False):
	
	"""Loops through all channels of image and returns channel of image with maximal range.

	Args:
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		tuple: Tuple containing:
			* img (numpy.ndarray): Monochromatic image.
			* ind_max (int): Index of channel chosen.
	"""
	
	#Empty list
	im_ranges=[]
		
	#Loop and save range	
	for k in range(len(np.shape(img))):
		im_ranges.append(img[:,:,k].max()-img[:,:,k].min())
	
	#Grab biggest range
	ind_max=im_ranges.index(max(im_ranges))
	
	#Reduce image
	img=img[:,:,ind_max]
	
	if debug:
		print "Warning, images are not monochromatic, going to choose channel number", ind_max+1,"!" 
	
	return img,ind_max

def loadImg(fn,enc,dtype='float'):
	
	"""Loads image from filename fn with encoding enc and returns it as with given dtype.

	Args:
		fn (str): File path.
		enc (str): Image encoding, e.g. 'uint16'.
			
	Keyword Args:
		dtype (str): Datatype of pixels of returned image.
	
	Returns:
		numpy.ndarray: Loaded image.
	"""
	
	#Load image
	img = skimage.io.imread(fn).astype(enc)
	
	#Getting img values
	img=img.real
	img=img.astype(dtype)
	
	return img

def normImg(img,imgPre,dataOffset=1.,debug=False):
	
	"""Norms image by preimage.

	Args:
		img (numpy.ndarray): Input image.
		imgPre (numpy.ndarray): Norming mask.
		
	Keyword Args:
		dataOffset (float): Offset used for norming.
		debug (bool): Print debugging messages and show debugging plots.
		
	Returns:
		numpy.ndarray: Processed image.
	"""
	
	#Norm and add offset to pre img to avoid singularities
	img=(img+dataOffset)/(imgPre+dataOffset)
			
	return img
	
def medianFilter(img,radius=1,debug=False,dtype='uint16',scaleImg=False,method='scipy',axes=None):
	
	"""Applies median filter to image.
	
	.. note:: Skimage algorithm requires images to be scaled up into 16bit range.
	
	Args:
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		radius (int): Odd integer defining median kernel size.
		axes (list): List of matplotlib axes used for plotting. If not specified, will generate new ones.
		debug (bool): Print debugging messages and show debugging plots.
		dtype (str): Optimal image dtype necessary for scaling.
		scaleImg (bool): Scale image to full range of dtype.
		method (str): Median algorithm used (scipy/skimage).
		
	Returns:
		numpy.ndarray: Processed image.
	"""
	
	if method=='skimage' and scaleImg==False:
		printWarning('Method skimage requires scaleImg=True. Going to do this.')  
		scaleImg=True
	
	if scaleImg:
		#Grab original dtype
		orgDtype=img.dtype
		
		#Grab original image
		orgImg=img.copy()
		
		#Grab maximum allowed value for dtype
		minval,maxval=getIntRangeDtype(dtype)
		
		#Scale up on range of uint 16 if necessary
		scale=maxval/img.max()
		img=scale*img
		
		#Convert to uint16
		img=img.astype(dtype)
		
	#Apply median rank filter
	if method not in ['scipy','skimage']:
		printWarning("Method " + method + " unknown. Going to use scipy(faster).")
	
	if method=='scipy':
		img=spsig.medfilt(img,kernel_size=int(radius))
	elif method=='skimage':
		img=skifilt.rank.median(img, skimage.morphology.disk(radius))
	
	if scaleImg:
		#Convert to old dtype	
		img=img.astype(orgDtype)
		
		#Scale back
		img=1./scale * img
		
	#Debugging plots
	if debug:
		
		#Make figure
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([2,2],titles=["Original Image", "After median","Histogram Original","Histogram median"],sup="medianFilter debugging output")	
		
		#Get common range
		vmin,vmax=getCommonRange([orgImg,img])
		
		showImgAndHist(orgImg,axes=axes[0:2],vmin=vmin,vmax=vmax)
		showImgAndHist(img,axes=axes[2:],vmin=vmin,vmax=vmax)
	
	return img
	
def gaussianFilter(img,sigma=2.,debug=False,axes=None):
	
	"""Applies gaussian filter to image.
	
	Args:
		img (numpy.ndarray): Input image.
			
	Keyword Args:
		sigma (float): Standard deviation of gaussian kernel applied.
		axes (list): List of matplotlib axes used for plotting. If not specified, will generate new ones.
		debug (bool): Print debugging messages and show debugging plots.
		
	Returns:
		numpy.ndarray: Processed image.
	"""
	
	
	#Grab original image
	orgImg=img.copy()
	
	#Apply gaussian filter
	img = skifilt.gaussian_filter(img,sigma)
	
	#Debugging plots
	if debug:
		
		#Make figure
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([2,2],titles=["Original Image", "After gaussian","Histogram Original","Histogram gaussian"],sup="gaussianFilter debugging output")	
		
		#Get common range
		vmin,vmax=getCommonRange([orgImg,img])
		
		showImgAndHist(orgImg,axes=axes[0:2],vmin=vmin,vmax=vmax)
		showImgAndHist(img,axes=axes[2:],vmin=vmin,vmax=vmax)
	
	return img

def substractBkgd(img,bkgd,substractMean=True,nonNeg=True):
	
	"""Substracts background from image.
	
	Args:
		img (numpy.ndarray): Input image.
		bkgd (numpy.ndarray): Background mask that will be used.
			
	Keyword Args:
		substractMean (bool): Substract mean background.
		nonNeg (bool): Fill negative pixels with 1. to avoid further problems.
	
	Returns:
		numpy.ndarray: Processed image.
	"""
	
	if substractMean:
		bkgd=np.mean(bkgd)
	
	#Make sure that there image stays non-negative
	if nonNeg:
		img,indX,indY=fixedThresh(img,0.,smaller=True,fill=1.)
	
	return img

def processImg(img,processDic,flatteningMask,bkgdMask,preMask,dataOffset=1.,axes=None,debug=False):

	"""Main image processing function containing the following steps:
	
	* Quadrant reduction.
	* Median filter.
	* Gaussian filter.
	* Norming.
	* Background substraction.
	* Flattening.
		
	Args:
		img (numpy.ndarray): Input image.
		processDic (dict): Dictionary defining what to do. See also pyfrp.subclasses.pyfrp_analysis .
		flatteningMask (np.ndarray): Flattening mask that will be used if flattening is selected.
		bkgdMask (numpy.ndarray): Background mask that will be used if background substraction is selected.
		preMask (numpy.ndarray): Preimage mask that will be used if norming is selected.
		
			
	Keyword Args:
		dataOffset (float): Offset used for norming.
		axes (list): List of matplotlib axes used for plotting. If not specified, will generate new ones.
		debug (bool): Print debugging messages and show debugging plots.
	
	Returns:
		numpy.ndarray: Processed image.
	"""
	
	
	#Flip image in case of quad_red and flip_before_process
	if 'quad' in processDic.keys():
		if 'flipBeforeProcess' in processDic.keys():
			img=flipQuad(img,debug=debug)
	
	#Apply median filter for denoising
	if 'median' in processDic.keys():
		img = medianFilter(img,radius=processDic['median'],debug=debug)
	
	#Apply gaussian blur to smooth out image
	if 'gaussian' in processDic.keys() and 'norm' not in processDic.keys():
		img = gaussianFilter(img,sigma=processDic['gaussian'],debug=debug)
	
	#Flatten img
	if 'flatten' in processDic.keys():
		img = flattenImg(img,flatteningMask)
	
	#Background Substraction
	if 'bkgd' in processDic.keys():
		img = substractBkgd(img,bkgdMask,substractMean=False)
	
	#Normalize by pre image
	if 'norm' in processDic.keys():
		img = normImg(img,preMask,dataOffset=dataOffset,debug=debug)
	
	#Quad reduction
	if 'quad' in processDic.keys():
		if 'flipBeforeProcess' in processDic.keys():
			#Flip image back in case of quad_red and flip_before_process (Need to unflip so that indices found by get_ind_regions still fit with image size)
			img=unflipQuad(img,debug=debug)
		else:	
			#Flip and unflip image again to average over first quadrant (Need to unflip so that indices found by get_ind_regions still fit with image size)
			img=flipQuad(img,debug=debug)
			img=unflipQuad(img,debug=debug)	
	
	if debug or axes!=None:
		
		#Build title
		
		###NOTE: need to change stuff here
		
		var=['quad_red','flip_before_process','norm_by_pre','gaussian']
		dic=pyfrp_misc_module.vars2dict(var,dict(locals()))
		title=pyfrp_misc_module.dict2string(dic,newline=True)
		
		if axes==None:
			
			#Create figure
			fig,axes = pyfrp_plot_module.makeSubplot([2,1],titles=title,sup="processImg debugging output")
		
		show_img_and_hist(img,axes=axes[0:2],title=[title,''])
		
	return img

def imgHist(img,binMin=0,binMax=65535,nbins=256,binSize=1,binsFit=True,fixSize=False,density=False):
	
	"""Creates histogram for image with some good default settings.
		
	Args:
		img (numpy.ndarray): Input image.
		
	Keyword Args:
		nbins (int): Number of bins of histogram.
		binSize (int): Size of bins.
		binMin (float): Minimum value of bins.
		binMax (float): Maximum value of bins.
		binsFit (bool): Fits bin range properly around range of image.
		fixSize (bool): Use binSize for fixed bin size.
		density (bool): Normalize histogram as probability density function.
		
	Returns:
		tuple: Tuple containing:

			* bins (float): Bin array of histogram.
			* hist (float): Histogram array.
			* w (float): Width of bins.
	"""
	
	if binsFit:
		if fixSize:
			bins=np.arange(0.9*np.nanmin(img),1.1*np.nanmax(img)+1,binSize)
		else:
			bins=np.linspace(np.nanmin(img)-1,np.nanmax(img)+1,nbins+1)
	else:
		bins=np.linspace(binMin,binMax,nbins+1)
	

	hist,binEdges=np.histogram(img,bins,density=density)
	

	bins=binEdges[:-1]+np.diff(binEdges)/2.
	
	w=(max(bins)-min(bins))/(2*len(bins))
	
	return bins,hist,w

def fixedThresh(img,thresh,smaller=False,fill=np.nan):
	
	"""Apply fixed threshhold to image and fill pixels with values 
	greater than thresh with fill value .
		
	Args:
		img (numpy.ndarray): Image for threshholding.
		thresh (float): Threshhold used.
		
	Keyword Args:
		fill (float): Fill values for pixels that are smaller or greater than thresh. 
		smaller (bool): Apply fill to pixels smaller or greater than thresh.
		
	Returns:
		tuple: Tuple containing:

			* img (numpy.ndarray): Output image.
			* indX (numpy.ndarray): x-indices of pixels that where threshholded. 
			* indY (numpy.ndarray): y-indices of pixels that where threshholded. 
	"""
	
	#Find pixels greater than thresh 
	if smaller:
		indX,indY=np.where(img<thresh)
	else:
		indX,indY=np.where(img>thresh)
	
	#fill with value
	img[indX,indY]=fill
	
	return img,indX,indY

def showImgAndHist(img,axes=None,sup='',title=None,color='b',vmin=None,vmax=None,binMin=0,binMax=65535,nbins=256,binSize=1,binsFit=True,fixSize=False,density=False):
	
	"""Creates plot of image and corresponding histogram.
		
	Args:
		img (numpy.ndarray): Image to be plotted.
		
	Keyword Args:
		nbins (int): Number of bins of histogram.
		binSize (int): Size of bins.
		binMin (float): Minimum value of bins.
		binMax (float): Maximum value of bins.
		binsFit (bool): Fits bin range properly around range of image.
		fixSize (bool): Use binSize for fixed bin size.
		density (bool): Normalize histogram as probability density function.
		vmin (float): Minimum display value of image.
		vmax (float): Maximum display value of image.
		title (list): List of titles of plots.
		sup (str): Supporting title of figure.
		axes (list): List of matplotlib axes used for plotting. If not specified, will generate new ones.
		color (str): Color of histogram plot.
		
	Returns:
		tuple: Tuple containing:

			* bins (float): Bin array of histogram.
			* hist (float): Histogram array.
			* w (float): Width of bins.
	"""
	
	if axes!=None:
		if len(axes)!=2:
			printWarning("Axes do not have right size! Will create new figure.")
			axes=None
		
	if title!=None:
		if len(title)!=2:
			printWarning("title do not have right size! Will not draw title.")
			title=None
	
	if axes==None:
		#Create figure
		fig,axes = pyfrp_plot_module.makeSubplot([2,1],titles=title,sup=sup,tight=True)
		
	#Show img
	a=axes[0].imshow(img,vmin=vmin,vmax=vmax)
	
	#Show histogram
	bins,hist,w=imgHist(img,binMin=binMin,binMax=binMax,nbins=nbins,binSize=binSize,binsFit=binsFit,fixSize=fixSize,density=density)
	axes[1].plot(bins,hist,color=color)
	
	if title!=None:
		if len(title)==2:
			axes[0].set_title(title[0])
			axes[1].set_title(title[1])	
	
	#Draw
	plt.draw()
	
	return bins,hist,w	
		
def getCommonRange(imgs):
	
	"""Finds common range of list of images.
		
	Args:
		imgs (list): List of images.
		
	Returns:
		tuple: Tuple containing:

			* minVal (float): Lowest value over all images.
			* maxVal (float): Largest value over all images.
	
	"""
	
	
	
	mins=[]
	maxs=[]
	
	for img in imgs:
		mins.append(img.min())
		maxs.append(img.max())
	
	minVal=min(mins)
	maxVal=max(mins)
	
	return minVal,maxVal

def getIntRangeDtype(dtype):
	
	"""Returns range of int based dtype.
	"""
	
	return np.iinfo(dtype).min,np.iinfo(dtype).max
	
def flattenImg(img,mask):
	
	"""Flattens image with flattening mask.
	"""
	
	return img*mask

def computeFlatMask(img,dataOffset):
	
	"""Computes flattening mask from image.
		
	Args:
		img (numpy.ndarray): Image to be used.
		dataOffset (int): Offset used.
			
	Returns:
		numpy.ndarray: Flattening mask.	
	"""
	
	mask=(img.max()+dataOffset)/(img+dataOffset)
	return mask

def computeMeanImg(fnFolder,fileList,dataEnc,median=False):
	
	"""Computes Mean Image from a list of files.
		
	Args:
		fnFolder (str): Path to folder containing files.
		fileList (list): List of file names in fnFolder.
		dataEnc (str): Encoding of images, e.g. uint16.
	
	Keyword Args:
		median (bool): Apply median filter per image.
		
	Returns:
		numpy.ndarray: Mean image.
	"""
	
	if len(fileList)==0:
		raise ValueError("There are no images in %s" %(fnFolder+fileList))
		
	for i in range(len(fileList)):
		
		fn=fnFolder+fileList[i]
		
		img=loadImg(fn,dataEnc)
		
		if median:
			img=medianFilter(img)
		
		if i>0:
			mImg=mImg+img
		else:
			mImg=img
	
	
	mImg=mImg/float(len(fileList))
	return mImg

def getMeanIntensitiesImgs(fnFolder,fileList,dataEnc):
	
	"""Reads all images in folder, returns mean intensity vector.
		
	Args:
		fnFolder (str): Path to folder containing files.
		fileList (list): List of file names in fnFolder.
		dataEnc (str): Encoding of images, e.g. uint16.
	
	Returns:
		list: Array of mean intensities per image.
	"""
	
	meanIntensities=[]
	
	for i in range(len(fileList)):
	
		fn=fnFolder+fileList[i]

		img=loadImg(fn,dataEnc)
		
		meanIntensities.append(np.mean(img))
		
	return meanIntensities	
				
def radialImgHist(img,center,nbins=10,byMean=True,maxR=None):
	
	"""Computes radial histogram of image from center.
		
	Args:
		img (numpy.ndarray): Image to be profiled
		center (list): Center of image.
	
	Keyword Args:
		nbins (int): Number of bins of histogram.
		byMean (bool): Norm bins by number of items in bin.
		maxR (float): Maximum radius considered.
		
	Returns:
		tuple: Tuple containing:

			* bins (numpy.ndarray): Bin vector.
			* binsMid (numpy.ndarray): Array of midpoints of bins.
			* histY (numpy.ndarray): Array of number of items per bin.
			* binY (numpy.ndarray): Array of average value per bin.
		
	"""
	
	#Get resolution
	res=img.shape[0]
	
	#Compute max distance from corners
	if maxR==None:
		maxR=max(dist(center,[0,0]),dist(center,[res,res]))
	
	#Create vector with bins
	bins=np.linspace(0,maxR,nbins+1)
	
	#Create result vectors
	histY=np.zeros(nbins)
	binY=np.zeros(nbins)
	
	for i in range(res):
		for j in range(res):
			for k in range(nbins):
				d=dist([i,j],center)
				
				if k==0:
					if bins[k]<=d and d<=bins[k+1]:
						binY[k]=binY[k]+img[i,j]
						histY[k]=histY[k]+1
						break
				else:
					if bins[k]<d and d<=bins[k+1]:
						binY[k]=binY[k]+img[i,j]
						histY[k]=histY[k]+1
						break
		
	if byMean:
		binY=binY/histY
	
	binsMid = [(bins[j+1] + bins[j])/2. for j in range(len(bins)-1)]
	
	return bins,binsMid,histY,binY


def plotRadialHist(img,center,nbins=10,byMean=True,axes=None,color='r',linestyle='-',fullOutput=False,maxR=None,plotBinSize=False,label='',legend=False,linewidth=1.):
	
	"""Plots radial histogram of image from center.
		
	Args:
		img (numpy.ndarray): Image to be profiled
		center (list): Center of image.
	
	Keyword Args:
		nbins (int): Number of bins of histogram.
		byMean (bool): Norm bins by number of items in bin.
		maxR (float): Maximum radius considered.
		plotBinSize (bool): Include number of items in bin on secondary axis.
		ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
		color (str): Color of plot.
		linestyle (str): Linestyle of plot.
		fullOutput (bool): Also return result arrays from radialImgHist ?
		linewidth (float): Linewidth of plot.
		
	
	Returns:
		tuple: Tuple containing:
		
			* ax (matplotlib.axes): Matplotlib axes.
			* ax2 (matplotlib.axes): Matplotlib axes for secondary axes. None if plotBinSize==False
			* bins (numpy.ndarray): Bin vector.
			* binsMid (numpy.ndarray): Array of midpoints of bins.
			* histY (numpy.ndarray): Array of number of items per bin.
			* binY (numpy.ndarray): Array of average value per bin.
		
	"""
	
	bins,binsMid,histY,binY=radialImgHist(img,center,nbins=nbins,byMean=byMean,maxR=maxR)
	
	if axes==None:
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Histogram"],sup="Histogram")
		ax=axes[0]
		if plotBinSize:
			ax2=ax.twinx()
		
	elif isinstance(axes,list): 
		if len(axes)==1:
			ax=axes[0]
			if plotBinSize:
				ax2=ax.twinx()
		else:
			ax=axes[0]
			ax2=axes[1]	
	else:
		ax=axes
		if plotBinSize:
			ax2=ax.twinx()
	
	if not plotBinSize:
		ax2=None
	
	ax.plot(binsMid,binY,color=color,linestyle=linestyle,label=label,linewidth=linewidth)
	
	if linestyle=='-':
		newstyle='--'
	else:
		newstyle='-'
	
	if plotBinSize:
		ax2.plot(binsMid,histY,color=color,linestyle=newstyle,linewidth=linewidth)
	
	plt.draw()
	
	if legend:
		plt.legend()
		
	if fullOutput:
		return ax,ax2,binsMid,histY,binY
	else:
		return ax,ax2

def plotRadialProfile(img,center,ax=None,color='r',linestyle='-',fullOutput=False,linewidth=1.):
	
	"""Plots radial profile of image from center.
		
	Args:
		img (numpy.ndarray): Image to be profiled
		center (list): Center of image.
	
	Keyword Args:
		ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
		color (str): Color of plot.
		linestyle (str): Linestyle of plot.
		fullOutput (bool): Also return result arrays from computeRadialProfile ?
		linewidth (float): Linewidth of plot.
		
	Returns:
		tuple: Tuple containing:
		
			* ax (matplotlib.axes): Matplotlib axes.
			* r (numpy.ndarray): Array of radii.
			* v (numpy.ndarray): Array of corresponding image values.
		
	"""
	
	r,v=computeRadialProfile(img,center)
	
	if ax==None:
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Profile"],sup="Profile")
		ax=axes[0]
		
	ax.plot(r,v,color=color,linewidth=linewidth)
	
	plt.draw()
	
	if fullOutput:
		return ax,r,v
	else:
		return ax
	
def computeRadialProfile(img,center):
	
	"""Computes radial profile of image from center.
		
	Args:
		img (numpy.ndarray): Image to be profiled
		center (list): Center of image.
			
	Returns:
		tuple: Tuple containing:
		
			* r (numpy.ndarray): Array of radii.
			* v (numpy.ndarray): Array of corresponding image values.
	"""
	
	#Empty lists
	r=[]
	v=[]
	
	#Get resolution
	res=img.shape[0]
	
	#Compute radial profile
	for i in range(res):
		for j in range(res):
			r.append(dist([i,j],center))
			v.append(img[i,j])
	
	#Sort according to radius
	s=sorted(zip(r, v), key=lambda r:r[0])
	s=np.asarray(s)
	r=s[:,0]		
	v=s[:,1]
	
	return r,v

	
def dist(p1,p2):
	
	"""Computes euclidean distance between two points p1 and p2.
	
	"""
	
	return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
	

		
def findProblematicNormingPixels(img,imgPre,dataOffset,axes=None,debug=False):
	
	"""Checks which pixels are problematic for norming.
	
	Note that function will temporarily set 
	>>> np.errstate(divide='raise')
	so numpy RuntimeWarning actually gets raised instead of just printed.
	
	Args:
		img (numpy.ndarray): Image to be normed.
		imgPre (numpy.ndarray): Image to be normed.
		dataOffset (int): Offset used for norming.
		
	Keyword Args:
		axes (list): List with matplotlib axes used for plotting. If not specified,
		will generate new ones.
		debug (bool): Show debugging outputs and plots.
		
	Returns:
		numpy.ndarray: Binary mask of problematic pixels.
	"""
	
	pxs=np.zeros(np.shape(img))
	
	#Loop through img and see which pixels actually generate errors when dividing
	with np.errstate(divide='raise'): #Need this line, so numpy RuntimeWarning actually gets raised instead of just printed
		for i in range(img.shape[0]):
			for j in range(img.shape[1]):
			
				try:
					a=(img[i,j]+dataOffset)/(imgPre[i,j]+dataOffset)
				except FloatingPointError:
					pxs[i,j]=1
					
	#Make some nice debugging plots
	if debug:
		vmin=min([img.min(),imgPre.min()])
		vmax=max([img.max(),imgPre.max()])
		
		
		if axes==None:
			fig,axes=pyfrp_plot_module.makeSubplot([2,4])
		if len(axes)<8:
			printWarning("Axes do not have the right size, will swiftly create some new ones for you.")
			fig,axes=pyfrp_plot_module.makeSubplot([2,4])
			
		showImgAndHist(img,axes=[axes[0],axes[4]],vmin=vmin,vmax=vmax)
		showImgAndHist(imgPre,axes=[axes[1],axes[5]],vmin=vmin,vmax=vmax)
		showImgAndHist(pxs*img,axes=[axes[2],axes[6]])
		showImgAndHist(pxs*imgPre,axes=[axes[3],axes[7]])
		
	return pxs	
		
def findMinOffset(fnFolder,fileList,dataEnc,oldOffset=None,defaultAdd=1.,debug=False):
	
	"""Simple function that loops through all images in file list and returns minimum integer 
	that needs to be added such that all pixels are positiv. 
	
	Args:
		fnFolder (str): Path to folder containing files.
		fileList (list): List of file names in fnFolder.
		dataEnc (str): Encoding of images, e.g. uint16.
		
	Keyword Args:
		oldOffset (int): Take some other offset into account.
		defaultAdd (int): Default value added to minimal offset.
		debug (bool): Show debugging outputs.
		
	Returns:
		int: Minimal Offset
	"""
	
	#Check if there are images
	if len(fileList)==0:
		printError("fileList contains no files. Going return min(["+str(oldOffset)+","+str(defaultAdd)+"]")
		if oldOffset!=None:
			return min([oldOffset,defaultAdd])
		else:
			return defaultAdd
	
	#Loop through images and get minima
	mins=[]
	
	for i in range(len(fileList)):
		
		#Load Img
		fn=fnFolder+fileList[i]
		img=loadImg(fn,dataEnc)	
		
		mins.append(img.min())
	minVal=min(mins)
	
	#Some debugging output
	if debug:
		print "Minimum value over all images in ", fnFolder ," = ", minVal
		print "Proposed offset is hence = ", abs(minVal)+defaultAdd
		
		if oldOffset!=None:
			print "Old offset was ", oldOffset
			if oldOffset>abs(minVal)+defaultAdd:
				print "Going to return ", oldOffset
			else:
				print "Going to return ", abs(minVal)+defaultAdd
	
	if oldOffset!=None:
		if oldOffset>abs(minVal)+defaultAdd:
			return oldOffset

	return abs(minVal)+defaultAdd

def otsuImageJ(img,maxVal,minVal,debug=False):
	
	"""Python implementation of Fiji's Otsu algorithm. 
	
	See also http://imagej.nih.gov/ij/source/ij/process/AutoThresholder.java.

	Args:
		img (numpy.ndarray): Image as 2D-array.
		maxVal (int): Value assigned to pixels above threshhold.
		minVal (int): Value assigned to pixels below threshhold.
		
	Keyword Args:
		debug (bool): Show debugging outputs and plots.
		
	Returns:
		tuple: Tuple containing:
		
			* kStar (int): Optimal threshhold
			* binImg (np.ndarray): Binary image
	"""
	
	#Initialize values
	#L = img.max()
	L = 256
	S = 0 
	N = 0
	
	#Compute histogram
	data,binEdges=np.histogram(img,bins=L)
	binWidth=np.diff(binEdges)[0]
	
	#Debugging plot for histogram
	if debug:
		binVec=arange(L)
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(121)
		ax.bar(binVec,data)
		plt.draw()
			
	for k in range(L):
		#Total histogram intensity
		S = S+ k * data[k]
		#Total number of data points
		N = N + data[k]		
	
	#Temporary variables
	Sk = 0
	BCV = 0
	BCVmax=0
	kStar = 0
	
	#The entry for zero intensity
	N1 = data[0] 
	
	#Look at each possible threshold value,
	#calculate the between-class variance, and decide if it's a max
	for k in range (1,L-1): 
		#No need to check endpoints k = 0 or k = L-1
		Sk = Sk + k * data[k]
		N1 = N1 + data[k]

		#The float casting here is to avoid compiler warning about loss of precision and
		#will prevent overflow in the case of large saturated images
		denom = float(N1 * (N - N1)) 

		if denom != 0:
			#Float here is to avoid loss of precision when dividing
			num = ( float(N1) / float(N) ) * S - Sk 
			BCV = (num * num) / denom
		
		else:
			BCV = 0

		if BCV >= BCVmax: 
			#Assign the best threshold found so far
			BCVmax = BCV
			kStar = k
	
	kStar=binEdges[0]+kStar*binWidth
	
	#Now manipulate the image
	binImg=np.zeros(np.shape(img))
	for i in range(np.shape(img)[0]):
		for j in range(np.shape(img)[1]):
			if img[i,j]<=kStar:
				binImg[i,j]=minVal
			else:				
				binImg[i,j]=maxVal
	
	if debug:
		print "Optimal threshold = ", kStar
		print "#Pixels above threshold = ", sum(binImg)/float(maxVal)
		print "#Pixels below threshold = ", np.shape(img)[0]**2-sum(binImg)/float(maxVal)
		
		ax2=fig.add_subplot(122)
		ax2.contourf(binImg)
		plt.draw()
		raw_input()
			
	return kStar,binImg	

def extractMicroscope(folder,ftype,fijiBin=None,macroPath=None):
	
	"""Converts all microscopy files of type ftype in folder to files using Fiji.

	Args:
		folder (str): Path to folder containing czi files
		ftype (str): Type of microscopy file, such as lsm or czi
	
	Keyword Args:
		fijiBin (str): Path to fiji binary
		macroPath (str): Path to fiji macro
		
	Returns:
		int: Returns 0 if success, -1 if error

	"""
	
	if ftype=='lsm':
		r=extractLSM(folder,fijiBin=fijiBin,macroPath=macroPath)
	elif ftype=='czi':
		r=extractCZI(folder,fijiBin=fijiBin,macroPath=macroPath)
	else:
		printError("Unknown filetype "+ftype)
		return -1
	return r

def extractLSM(folder,fijiBin=None,macroPath=None):
	
	"""Converts all lsm files in folder to tif files using Fiji.

	Args:
		folder (str): Path to folder containing lsm files
	
	Keyword Args:
		fijiBin (str): Path to fiji binary
		macroPath (str): Path to fiji macro
		
	Returns:
		int: Returns 0 if success, -1 if error

	"""

	if macroPath==None:
		macroPath=pyfrp_misc_module.getMacroDir()+'lsm2tif.ijm'
	
	return runFijiMacro(macroPath,folder,fijiBin=fijiBin)
	
def extractCZI(folder,fijiBin=None,macroPath=None):
	
	"""Converts all czi files in folder to tif files using Fiji.

	Args:
		folder (str): Path to folder containing czi files
		
	Keyword Args:	
		fijiBin (str): Path to fiji binary
		macroPath (str): Path to fiji macro
		
	Returns:
		int: Returns 0 if success, -1 if error

	"""
	
	if macroPath==None:
		macroPath=pyfrp_misc_module.getMacroDir()+'czi2tif.ijm'
	
	return runFijiMacro(macroPath,folder,fijiBin=fijiBin)
	
	
def runFijiMacro(macroPath,macroArgs,fijiBin=None):
	
	"""Runs Fiji Macro.

	Args:
		macroPath (str): Path to fiji macro
		macroArgs (str): Arguments being passed to Fiji macro
		
	Keyword Args:	
		fijiBin (str): Path to fiji binary
			
	Returns:
		int: Returns 0 if success, -1 if error

	"""
	
	if fijiBin==None:
		fijiBin=pyfrp_misc_module.getFijiBin()
		
	#Define Command to run
	cmd=fijiBin+" -macro "+ macroPath + " '"+macroArgs +"'"+ " -batch " 
	
	#Run
	try:
		os.system(cmd)
		return 0
	except:
		printError("Something went wrong executing:" )
		print cmd
		return -1	
	
	
	
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Functions that still need testing
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Show all processing possibilities for image (Flip/Norm/Gauss)

#def show_process_poss(img,fn_pre,dataEnc):
	
	##Create figure and axes
	#fig,axes1 = pyfrp_plot_module.makeSubplot([2,4],tight=False)
	
	##Run all possibilities with quad_red=False
	#img1=processImg(img,fn_pre,dataEnc,quad_red=False,flip_before_process=True,norm_by_pre=False,gaussian=False,gaussianSigma=2.,dataOffset=1.,axes=[axes1[0],axes1[4]],debug=False)
	#img2=processImg(img,fn_pre,dataEnc,quad_red=False,flip_before_process=True,norm_by_pre=True,gaussian=False,gaussianSigma=2.,dataOffset=1.,axes=[axes1[1],axes1[5]],debug=False)
	#img3=processImg(img,fn_pre,dataEnc,quad_red=False,flip_before_process=True,norm_by_pre=False,gaussian=True,gaussianSigma=2.,dataOffset=1.,axes=[axes1[2],axes1[6]],debug=False)
	#img4=processImg(img,fn_pre,dataEnc,quad_red=False,flip_before_process=True,norm_by_pre=True,gaussian=True,gaussianSigma=2.,dataOffset=1.,axes=[axes1[3],axes1[7]],debug=False)
	
	##Create figure and axes
	#fig,axes2 = pyfrp_plot_module.makeSubplot([2,4],tight=False)
	
	##Run all possibilities with quad_red=True and flip_before_process=True
	#img5=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=True,norm_by_pre=False,gaussian=False,gaussianSigma=2.,dataOffset=1.,axes=[axes2[0],axes2[4]],debug=False)
	#img6=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=True,norm_by_pre=True,gaussian=False,gaussianSigma=2.,dataOffset=1.,axes=[axes2[1],axes2[5]],debug=False)
	#img7=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=True,norm_by_pre=False,gaussian=True,gaussianSigma=2.,dataOffset=1.,axes=[axes2[2],axes2[6]],debug=False)
	#img8=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=True,norm_by_pre=True,gaussian=True,gaussianSigma=2.,dataOffset=1.,axes=[axes2[3],axes2[7]],debug=False)
	
	##Create figure and axes
	#fig,axes3 = pyfrp_plot_module.makeSubplot([2,4],tight=False)
	
	##Run all possibilities with quad_red=True and flip_before_process=False
	#img9=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=False,norm_by_pre=False,gaussian=False,gaussianSigma=2.,dataOffset=1.,axes=[axes3[0],axes3[4]],debug=False)
	#img10=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=False,norm_by_pre=True,gaussian=False,gaussianSigma=2.,dataOffset=1.,axes=[axes3[1],axes3[5]],debug=False)
	#img11=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=False,norm_by_pre=False,gaussian=True,gaussianSigma=2.,dataOffset=1.,axes=[axes3[2],axes3[6]],debug=False)
	#img12=processImg(img,fn_pre,dataEnc,quad_red=True,flip_before_process=False,norm_by_pre=True,gaussian=True,gaussianSigma=2.,dataOffset=1.,axes=[axes3[3],axes3[7]],debug=False)
	
	##Adjust display ranges
	#imgs_norm=[img2,img4,img6,img8,img10,img12]
	#imgs=[img1,img3,img5,img7,img9,img11]
	
	#axes_norm=[axes1[1],axes1[3],axes2[1],axes2[3],axes3[1],axes3[3]]
	#axes_org=[axes1[0],axes1[2],axes2[0],axes2[2],axes3[0],axes3[2]]
	
	#min_norm,max_norm=get_common_range(imgs_norm)
	#min_org,max_org=get_common_range(imgs)
	
	#pyfrp_plot_module.adjust_imshow_range(axes_org,vmin=min_org,vmax=max_org)
	#pyfrp_plot_module.adjust_imshow_range(axes_norm,vmin=min_norm,vmax=max_norm)
	
	#print "done"
	#raw_input()



##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Find longest boundary 

#def findLongestBoundary(vals,max_val):
	
	##Find longest boundary contour
	#max_lengths=[]
	
	#for i in range(max_val):
		#conts_embr=skimsr.find_contours(vals, i, fully_connected='high', positive_orientation='low')
		
		#if np.shape(conts_embr)[0]>0:
		
			#lengths=[]
			#for cont in conts_embr:
				#lengths.append(len(cont))
			
			#max_lengths.append(max(lengths))
		
		#else:
			#max_lengths.append(-100)
		
	#max_len=max(max_lengths)
	#opt_thresh=max_lengths.index(max_len)
	
	##Grab again longest boundary
	#conts_final=skimsr.find_contours(vals, opt_thresh, fully_connected='high', positive_orientation='low')
	
	##Extracting longest boundary
	#cont_lengths=[]
	#for i in range(np.shape(conts_final)[0]):
		#cont_lengths.append(np.shape(conts_final[i])[0])
	
	#indvec=find_range_ind(cont_lengths,20,max_len+10)
	
	#conts_final_new=[]
	#for i in indvec:
		#conts_final_new.append(conts_final[i])	
	
	#return conts_final_new, max_len, opt_thresh

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Compute center of mass

#def computeCOM(vals):
	
	##Compute center of mass
	#mass=0
	#x_weights=0
	#y_weights=0
	
	#for i in range(np.shape(vals)[0]):
		#for j in range(np.shape(vals)[1]):
			#mass=mass+vals[i,j]
			#x_weights=x_weights+i*vals[i,j]
			#y_weights=y_weights+j*vals[i,j]
			
	#com=[x_weights/mass,y_weights/mass]		
	
	#return com

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Compute final radius through histogram

#def computeRadiusHist(radiuses,debug_opt):
	
	#bins=round(len(radiuses)/20)
	#data,bin_edges=histogram(radiuses,bins=bins)
	
	#bin_vec=linspace(min(radiuses),max(radiuses),bins)
	
	#if debug_opt==1:
		#fig=plt.figure()
		#fig.show()
		#ax=fig.add_subplot(111)
		#ax.bar(bin_vec,data)
		
		#plt.draw()
		#raw_input()
		
	#ind=data.argmax()
	
	#return radiuses[ind]
		
##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##find index of all entries in vec with vals between valmin and valmax

#def find_range_ind(vec,valmin,valmax):
	#indvec=[]
	#for i in range(np.shape(vec)[0]):
		#if vec[i] >= valmin and vec[i] <= valmax:
			#indvec.append(i)
			
	#return indvec	
	
##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Crop concentration of rim from pre image

#def cropConcRimFromPre(embryo):
	
	##Open image file
	#fnImg=embryo.fn_preimage
	#fnImg=str(fnImg)
	
	##Reading in image
	#im = skimage.io.imread(fnImg).astype(embryo.dataEnc)
	
	##Getting img values
	#imgVals=im.real
	
	##Get indices of domains
	#ind_circ_x,ind_circ_y,ind_sq_x,ind_sq_y,ind_slice_x,ind_slice_y, ind_rim_x,ind_rim_y=get_ind_regions(embryo.offset_bleached_px,embryo.side_length_bleached_px,embryo.radius_embr_px,embryo.center_embr_px,embryo.rim,embryo.add_rim_from_radius,imgVals,0)
	
	#num_in_slice=np.shape(ind_slice_x)[0]
	
	##Calculate concentration in slice
	#c_slice,c_rim=analyze_conc_slice(ind_slice_x,ind_slice_y,imgVals,n0,ind_rim_x,ind_rim_y,0,embryo.add_rim_img)
	
	##Save to embryo object
	#embryo.conc_rim=c_slice
	
	#return embryo

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Smooth out image via sliding window

#def slidingWin(embryo):
	
	##Check if overlap options is turned on
	#if embryo.slide_overlap==0:
		
		##Check if slide_win_width is a divisor of data_res_px, if not, adjust it
		#if mod(embryo.data_res_px,embryo.slide_win_width)!=0:
			#if embryo.debug_analysis==1:
				#print "Sliding window size is not divisor of img size, going to correct that"
			
			#d=float("inf")
			
			##Getting all divisors of data_res_px
			#for i in range(1,int(embryo.data_res_px/2+1)):
				
				##Check if divisor
				#if mod(embryo.data_res_px,i)==0:
					
					##Compute distance to slide_win_width
					#d_new=abs(embryo.slide_win_width-i)
					#if d_new<d:
						##Divisor that lies closer to slide_win_width
						#i_new=i
						#d=d_new
					#else:
						##Won't get better, so break
						#break
			
			#if embryo.debug_analysis==1:
				#print "Old slide_win_width:", embryo.slide_win_width, "New slide_win_width", i_new
			
			##Set new slide_win_width
			#embryo.slide_win_width=i_new
		
		##Calculate number of slide_steps
		#slide_steps=int(embryo.data_res_px/embryo.slide_win_width)
		
	#elif embryo.slide_overlap==1:
		##Calculate number of slide_steps
		#slide_steps=int(embryo.data_res_px-(embryo.slide_win_width-1))
	
	##Slide matrix
	#slides=zeros((slide_steps,slide_steps))
	
	#if embryo.slide_img=="post":
		#slide_img=embryo.im_reg_ICs
	#elif embryo.slide_img=="pre":	
		#fnImg=embryo.fn_preimage
		#fnImg=str(fnImg)
		#im = skimage.io.imread(fnImg).astype(embryo.dataEnc)
		
		#slide_img=im.real
		#slide_img=slide_img.astype('float')
	
	
	##Go through all slides and crop current slide from img and calculate mean
	#for i in range(slide_steps):
		#for j in range(slide_steps):
			#if embryo.slide_overlap==0:
				#curr_win=slide_img[embryo.slide_win_width*i:embryo.slide_win_width*(i+1),embryo.slide_win_width*j:embryo.slide_win_width*(j+1)]
			#elif embryo.slide_overlap==1:
				#curr_win=slide_img[i:i+embryo.slide_win_width,j:embryo.slide_win_width+j]
				
			#slides[i,j]=mean(curr_win)
			
	##Some debugging plots
	#if embryo.debug_analysis==1:
	
		#fig=plt.figure()
		#fig.show()
		#cflevels=linspace(slide_img.min(),slide_img.max(),10)
		#ax=fig.add_subplot(121)
		#ax.contourf(slide_img,levels=cflevels)
		#ax=fig.add_subplot(122)
		#ax.contourf(slides,levels=cflevels)
		
		#plt.draw()
		#raw_input()
	
	#embryo.slides=slides
			
	#return embryo			