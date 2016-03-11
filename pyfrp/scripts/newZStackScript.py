#Script to analyze zstack images and compute geomerty

#=========================================================================================================================
#Importing modules
#=========================================================================================================================

import pyfrp_img_module
import pyfrp_misc_module
import pyfrp_gmsh_module
import pyfrp_zstack_module


import scipy.interpolate as interp



#import skimage.io as skiio
#import skimage.measure as skimsr
import matplotlib.pyplot as plt
#import cv2
#from scipy import ndimage
#import os
import numpy as np

#=========================================================================================================================
#Parameters
#=========================================================================================================================

#~~~~~~~~~~~~~~~~~~
#Files

#Get current zstack
folder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/20160119_zstack_nobeads/01"
files=pyfrp_misc_module.getSortedFileList(folder,".tif")
outfile="zStack2"
outpath="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/outfiles/"

#Some flags what to do
debug=True
debugFill=True
sumPlot=True
threshMeth='otsu'

reverseOrder=True
composeContours=False

#Some empty vectors to save results
mpolys=[]

#Some parameters
stackDistance=4.51
volSizePx=25
fillSamples=2000
fillMode="random"
dPoly=1

#Reconstruction algorithm
reconAlg="poisson"
#recon_alg="vcg"


meshlabscript="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/scripts/normals_"+reconAlg+".mlx"

#=========================================================================================================================
#Automatically Detecting Geometry Boundaries
#=========================================================================================================================

#Creating figure for final plot
if sumPlot==1:
	figSum=plt.figure()
	figSum.show()
numRows=np.ceil(np.floor(len(files))/4.)
j=1

#If reverse_order is selected, reverse order of files
if reverseOrder:
	files.reverse()

#files=files[5:]

#Loop through files
for j,fn in enumerate(files):
	
	print "---------------------------------"
	print "Img:", fn
	
	#Reading img file
	img=pyfrp_img_module.loadImg(folder+"/"+fn,'uint16')
	
	#Apply gaussian filter
	img=pyfrp_img_module.gaussianFilter(img,sigma=2)
	
	#Applying otsu algorithm
	if threshMeth=='otsu':
		threshhold,thresh=pyfrp_img_module.otsuImageJ(img,1,0,0)
	elif threshMeth=='cutOff':
		threshhold=0.3*img.max()
		thresh=img
		thresh[(where(img<threshhold))]=0		
	else:
		threshhold=0.76*img.max()
		thresh=img
	
	#Find contours
	contours,thresh,areas=pyfrp_zstack_module.getContours(thresh)
	
	#Grab biggest area contour
	idxMax=areas.index(max(areas))
	contourMax=contours[idxMax]
	
	#To np array
	contours=np.asarray(contours)
	
	#Grab contours
	contoursFinal=[]
	if len(np.shape(contours))==4:
		contoursFinal.append(contours[0,:,0,:])
	else:
		for c in contours:
			contoursFinal.append(c[:,0,:])
	
	#Plot
	if sumPlot:
		axSum=figSum.add_subplot(4,numRows,j)
		axSum.imshow(img,cmap="Greys")
		axSum.set_title(fn)
		for c in contoursFinal:
			axSum.plot(c[:,0],c[:,1],'r')
		
		contourMax=np.asarray(contourMax)
		axSum.plot(contourMax[:,0,0],contourMax[:,0,1],'b')
		plt.draw()
		
