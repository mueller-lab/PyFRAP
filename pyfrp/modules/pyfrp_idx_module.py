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

#Indexing module for PyFRAP toolbox, including following functions:

#(1)  getCircleIdxImg 
#(2)  getRectangleIdxImg
#(3)  getSquareIdxImg
#(4)  getAllIdxImg
#(5)  getPolyIdxImg
#(6)  getCircleIdxMesh
#(7)  getSliceIdxMesh
#(8)  getRectangleIdxImg
#(9)  getSquareIdxImg
#(10)  getPolyIdxMesh
#(11)  checkInsideCircle
#(12)  checkInsideSquare
#(13)  checkInsideRectangle
#(14)  checkInsidePoly
#(15)  checkQuad
#(16)  checkSquareSize
#(17)  checkSquareCenteredFromInd
#(18)  ind2quad
#(19)  regions2quad
#(20)  ind2mask

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
import os, os.path
import sys

#PyFRAP modules
import pyfrp_misc_module as pyfrp_misc
import pyfrp_plot_module as pyfrp_plt


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in image that lie within given circle

def getCircleIdxImg(center,radius,res,debug=False):
	
	#Empty index vectors
	ind_circ_x=[]
	ind_circ_y=[]
	
	if debug:
		ind_slice_debug=np.zeros((res,res))
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in circle
			if checkInsideCircle(i+1,j+1,center,radius):
				ind_circ_x.append(i)
				ind_circ_y.append(j)
				if debug:
					ind_circ_debug[i,j]=1
	
	if debug:
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Circle"],sup="getCircleIdxImg debugging output")
		axes[0].contourf(ind_circ_debug)
		
	return ind_circ_x,ind_circ_y

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in image that lie within given rectangle

def getRectangleIdxImg(offset,sidelengthX,sidelengthY,res,debug=False):
	
	indX=[]
	indY=[]
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in square
			if checkInsideRectangle(i+1,j+1,offset,sidelengthX,sidelengthY):
				indX.append(i)
				indY.append(j)
				
				if debug:
					indDebug[i,j]=1
					
	if debug:
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Rectangle"],sup="getRectangleIdxImg debugging output")
		axes[0].contourf(ind_circ_debug)			
		
	return indX,indY	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in image that lie within given square

def getSquareIdxImg(offset,sidelength,res,debug=False):
	
	indX=[]
	indY=[]
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in square
			if checkInsideSquare(i+1,j+1,offset,sidelength):
				indX.append(i)
				indY.append(j)
				
				if debug:
					indDebug[i,j]=1
					
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Square"],sup="getSquareIdxImg debugging output")
		axes[0].contourf(ind_circ_debug)			
		
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in image

def getAllIdxImg(res,debug=False):
	
	#Empty index vectors
	indX=[]
	indY=[]
	
	if debug:
		indDebug=zeros((res,res))
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in circle
			indX.append(i)
			indY.append(j)
			if debug:
				indDebug[i,j]=1
				
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["slice"],sup="getAllIdxImg debugging output")
		axes[0].contourf(indDebug)
		
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in image that lie within given polygon

def getPolyIdxImg(corners,res,debug=False):
	
	#Convert to np array if necessary
	corners=np.asarray(corners)
	
	#Define polygonial patch
	poly = ptc.Polygon(corners,edgecolor='r',facecolor=(1,0,0,.2),)
	
	#Create grid
	x_int=np.arange(1,res+1,1)
	y_int=np.arange(1,res+1,1)
	g = np.meshgrid(x_int, y_int)
	
	#Zip them into coordinate tuples
	coords = list(zip(*(c.flat for c in g)))
	
	#Check which point is inside
	pts = np.vstack([p for p in coords if poly.contains_point(p, radius=0)])

	indX,indY= pts[0,:],pts[1,:]
	
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in mesh that lie within given circle
		
def getCircleIdxMesh(center,radius,mesh,zmin="-inf",zmax="inf",debug=False):
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.cellCenters
	
	#Convert into complex numbers
	c=np.array([np.complex(xc,yc) for xc,yc in zip(x,y)])
	centerC=np.complex(center[0],center[1])
	
	#Get indices in Circle
	indCircle=np.where(np.abs(c-centerC)<radius)[0]

	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	

	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indCircle)
	
	return indFinal

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in mesh that lie within given slice
	
def getSliceIdxMesh(z,zmin,zmax,debug=False):
	

	
	indSlice=np.where((z<zmax) & (z > zmin))[0]
	return indSlice

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in mesh that lie within given rectangle

def getRectangleIdxMesh(sidelengthX,sidelengthY,offset,mesh,zmin="-inf",zmax="inf",debug=False):
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.cellCenters
	
	#Getting indices
	indSquare=np.where((offset[0]<x) & (x<offset[0]+sidelengthX) & (offset[1]<y) & (y<offset[1]+sidelengthY))[0]
	
	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	
	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indFinal)
	
	return indFinal

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in mesh that lie within given square

def getSquareIdxMesh(sidelength,offset,mesh,zmin="-inf",zmax="inf",debug=False):
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.cellCenters
	
	#Getting indices
	indSquare=np.where((offset[0]<x) & (x<offset[0]+sidelength) & (offset[1]<y) & (y<offset[1]+sidelength))[0]
	
	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	
	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indSquare)
	
	return indFinal

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns all indices in mesh that lie within given polygon

def getPolyIdxMesh(corners,mesh,zmin="-inf",zmax="inf",debug=False):
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.cellCenters
	
	#Bookkeeping list
	indPoly=[]
	
	#Loop through coordinates and check if inside
	for i in range(len(x)):
		if checkInsidePoly(x[i],y[i],corners):
			indPoly.append(i)
	
	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	
	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indPoly)
		
	return indFinal	
		

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is in circle with given radius and center

def checkInsideCircle(x,y,center,radius):
	return np.sqrt((x-center[0])**2+(y-center[1])**2)<radius
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is in square with given offset and sidelength
	
def checkInsideSquare(x,y,offset,sidelength):
	return x<=offset[0]+sidelength and offset[0]<=x and y<=offset[1]+sidelength and offset[1]<=y
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is in square with given offset and sidelength
	
def checkInsideRectangle(x,y,offset,sidelengthX,sidelengthY):
	return x<=offset[0]+sidelengthX and offset[0]<=x and y<=offset[1]+sidelengthY and offset[1]<=y	
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is in square with given offset and sidelength

def checkInsideImg(x,y,res,offset=[0,0]):
	return (x>offset[0]) * (x<offset[0]+res)*(y>offset[1]) * (y<offset[1]+res)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is in polyogn, checks first if vector or just value

def checkInsidePolyVec(x,y,poly):	
	try:
		len(x)
	except TypeError:
		return checkInsidePoly(x,y,poly)
	
	vec=[]
	for i in range(len(x)):
		vec.append(checkInsidePoly(x[i],y[i],poly))
	     
	return vec

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is in polyogn

def checkInsidePoly(x,y,poly):

	#Taken from http://www.ariel.com.au/a/python-point-int-poly.html
	n = len(poly)
	inside =False
		
	
	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y
	
	return inside

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if coordinate (x,y) is first quadrant

def checkQuad(x,y,res):
	return (res/2.-1<=x) & (res/2.-1<=y)
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if square has right size, returns True if so.

def checkSquareSize(ind_sq_x,ind_sq_y,sidelength):

	if floor(sidelength)**2==len(ind_sq_x) and floor(sidelength)**2==len(ind_sq_x) and floor(sidelength)**2==len(unique(ind_sq_x))*len(unique(ind_sq_y)):
		return True
	else:
		return False
		
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if square is centered in image from indices (first bracket should be zero and second too)
	
def checkSquareCenteredFromInd(ind_sq_x,ind_sq_y,res):
	return not bool((min(ind_sq_x)-(res-1-max(ind_sq_x)))+(min(ind_sq_y)-(res-1-max(ind_sq_y))))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if square is centered in image (need a correction by .5 because there is a difference between pixels and coordinates)
	
def checkSquCentered(offset,sidelength,res):
	return not bool(((res-sidelength)-2*(offset[0]-0.5)) and ((res-sidelength)-2*(offset[1]-0.5)))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Reduces indices found for whole domain to first quadrant

def idx2QuadImg(indX,indY,res,debug=False):
	
	#Convert to np array
	indX=np.asarray(indX)
	indY=np.asarray(indY)
	
	#Find indices of first quadrant
	inds=where((res/2.-1<=indX) & (res/2.-1<=indY))[0]
	
	#Assign new indices
	indXQuad=indX[inds]
	indYQuad=indY[inds]
	
	#Debugging plot if necessary
	if debug:
	
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,3],titles=["Original indices","Flipped indices","Original and flipped indices"],sup="ind2quad debugging output")
		
		#Plot original and flipped indices
		axes[0].plot(indX,indY,'y*')
		axes[1].plot(indXQuad,indYQuad,'r*')
		axes[2].plot(indX,indY,'y*')
		axes[2].plot(indXQuad,indYQuad,'r*')
		
		#Draw
		plt.draw()
	
	return indXQuad, indYQuad

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Reduces indices of regions specified in inds found for whole domain to first quadrant

def regions2quad(inds,res,debug=False):
	inds_quad=[]
	for ind in inds:
		ind_quad=ind2quad(ind[0],ind[1],res,debug=debug)
		inds_quad.append(ind_quad)	
	
	return inds_quad


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fills img at given positions with value  

def ind2mask(vals,ind_x,ind_y,val):
	
	vals[ind_x,ind_y]=val
	
	return vals

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Converts mask into indices

def mask2ind(mask,res):
	
	#To bool
	mask=mask.astype(bool)

	#idx grid
	x=np.arange(res)
	y=np.arange(res)	
	X,Y=np.meshgrid(x,y)

	#Slice idx grid
	idxX_new=X[mask].flatten().astype(int)
	idxY_new=Y[mask].flatten().astype(int)
	
	return idxX_new, idxY_new

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds theoretical pixels that could be filled up with rim concentration for a square

def getExtendedPixelsSquare(offset,sidelength,res,debug=False):
	
	indX=[]
	indY=[]
	
	x=np.arange(np.floor(offset[0]),np.ceil(offset[0]+sidelength))
	y=np.arange(np.floor(offset[1]),np.ceil(offset[1]+sidelength))
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds theoretical pixels that could be filled up with rim concentration for a rectangle

def getExtendedPixelsRectangle(offset,sidelengthX,sidelengthY,res,debug=False):
	
	x=np.arange(np.floor(offset[1]),np.ceil(offset[0]+sidelengthX))
	y=np.arange(np.floor(offset[1]),np.ceil(offset[1]+sidelengthY))
	
	indX=[]
	indY=[]
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res) and checkInsideRectangle(i,j,offset,sidelengthX,sidelengthY):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds theoretical pixels that could be filled up with rim concentration for a circle

def getExtendedPixelsCircle(center,radius,res,debug=False):
	
	x=np.arange(np.ceil(center[0]-radius),np.floor(center[0]+radius))
	y=np.arange(np.ceil(center[1]-radius),np.floor(center[1]+radius))
	
	indX=[]
	indY=[]
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res) and checkInsideCircle(i,j,center,radius):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds theoretical pixels that could be filled up with rim concentration for a polygon

def getExtendedPixelsPolygon(corners,res,debug=False):
	
	cornersNP=np.array(corners)
	
	xmax=cornersNP[:,0].max()
	xmin=cornersNP[:,0].min()
	ymax=cornersNP[:,1].max()
	ymin=cornersNP[:,1].min()
	
	x=np.arange(np.floor(xmin),np.ceil(xmax))
	y=np.arange(np.floor(ymin),np.ceil(ymax))
	
	
	indX=[]
	indY=[]
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res) and checkPolygon(i,j,corners):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds theoretical pixels that could be filled up with rim concentration for a list of ROIs

def getCommonExtendedPixels(ROIs,res,debug=False,procedures=None):
	
	xExtend,yExtend=getCommonXYExtend(ROIs,debug=debug)
	
	x=np.arange(np.floor(xExtend[0]),np.ceil(xExtend[1]))
	y=np.arange(np.floor(yExtend[0]),np.ceil(yExtend[1]))
	
	indX=[]
	indY=[]
	
	if procedures==None:
		procedures=np.ones(np.shape(ROIs))
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res):
				b=True
				for k,r in enumerate(ROIs):
					if 1+procedures[k]:
						b=b and r.checkXYInside(i,j)
					else:					
						b=b and not r.checkXYInside(i,j)
				if b:
					indX.append(i)
					indY.append(j)
	
	return indX,indY

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gets common x-y-extend for a list of ROIs

def getCommonXYExtend(ROIs,debug=False):
	xExtends,yExtends=[],[]
	for r in ROIs:
		xExtend,yExtend=r.computeXYExtend()
		xExtends.append(xExtend)
		yExtends.append(yExtend)
	
	xExtends=np.array(xExtends)
	yExtends=np.array(yExtends)
	
	xExtend=[xExtends[:,0].min(),xExtends[:,1].max()]
	yExtend=[yExtends[:,0].min(),yExtends[:,1].max()]
	
	if debug:
		print "xExtend = ", xExtend
		print "yExtend = ", yExtend
	
	return xExtend, yExtend

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Remove repeated indices tupels from index lists for images

def remRepeatedImgIdxs(idxX,idxY,debug=False):
	idx=zip(idxX,idxY)
	idx=pyfrp_misc.remRepeatsList(idx)
	idxX,idxY=pyfrp_misc.unzipLists(idx)
	return idxX,idxY

	
				
				
	
	
	