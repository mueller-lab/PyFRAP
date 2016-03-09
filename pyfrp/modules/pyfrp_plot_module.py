#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Plotting module for image analysis and simulation for PyFRAP toolbox, including following functions:
#(1)  geom_wire_plot: Produces wireframe plot of embryo geometry.
#(2)  surf_plot: Produces surface plot of given slice values.
#(3)  cont_plot: Produces contour plot of given slice values.
#(4)  conc_plot: Plots concentration profiles until step=step.
#(5)  all_plot: Produces figure with surf_plot, cont_plot and conc_plot as subfigures
#(6)  plot3Dmesh: Plot mesh in given region. Can also highlight given regions of the mesh. (slow)
#(7)  plot3Dtet: Plots a list of given tetrahedrons from a list of mesh nodes.
#(8)  plot3Dnodes: Plots nodes in given region similarly to plot3Dmesh.
#(9)  gen_masks: Generates masks of slice, square and out for plotting. Fills values out of region with fill_value (preferably 0 or NaN).

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#numpy
import numpy as np

#Plotting
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Misc
import sys
import os

#PyFRAP
import pyfrp_img_module
from pyfrp_term_module import *

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================


#===================================================================================================================================
#Opens figure with first post image and give click selector for boundaries
#===================================================================================================================================

class FRAPBoundarySelector():
	
	def __init__(self,embryo):
		
		#Passing embryo to class
		self.embryo=embryo	
		
		#Plot resolution
		self.dpi = 100
	
		#Some bookkeeping variables
		self.radiusPt=[]
		self.centerPt=[]
		self.centerMarker=[]
		
		#Results
		self.radius=None
		self.center=None
		
		#Creating figure and canvas
		self.createCanvas()
			
		#Plot image if existent
		self.showFirstDataImg()
		plt.show()
	
	def createCanvas(self):
			
		h=500/self.dpi
		v=500/self.dpi
		
		self.fig = plt.figure()
		self.fig.show()
		self.fig.set_size_inches(h,v,forward=True)
		
		self.canvas=self.fig.canvas
		
		self.ax = self.fig.add_subplot(111)
		self.ax.set_xlim([0, self.embryo.dataResPx])
		self.ax.set_ylim([0, self.embryo.dataResPx])
		
		#Connect to mouseclick event
		self.fig.canvas.mpl_connect('close_event', self.closeFigure)
		self.canvas.mpl_connect('button_press_event', self.getMouseInput)
		
		self.canvas.draw()
		
		return 
	
	def closeFigure(self,event):
		return self.center,self.radius
		
	def getEmbryo(self):	
		return self.embryo
	
	def getResults(self):
		return self.center,self.radius
	
	def drawCenterMarker(self):
		
		centerImg=[self.embryo.dataResPx/2.,self.embryo.dataResPx/2.]
		
		if len(self.centerMarker)>0:
			self.clearCenterMarker()
		else:	
			pt=ptc.Circle(centerImg,radius=3,fill=True,color='y')
			self.centerMarker.append(self.ax.add_patch(pt))
			
		self.fig.canvas.draw()
		
		return self.centerMarker
		
	def clearCenterMarker(self):
	
		for pt in self.centerMarker:
			pt.remove()
			
		self.fig.canvas.draw()	
			
		self.centerMarker=[]
	
	def drawCenter(self):
		
		if len(self.centerPt)>0:
			self.clearCenter()
	
		pt=ptc.Circle(self.center,radius=3,fill=True,color='r')
		self.centerPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.centerPt
	
	def clearCenter(self):
	
		for pt in self.centerPt:
			pt.remove()
			
		self.centerPt=[]
		
		self.fig.canvas.draw()
		
	def drawRadius(self):
		
		if len(self.radiusPt)>0:
			self.clearRadius()
			
		pt=ptc.Circle(self.center,radius=self.radius,fill=False,color='r')
		self.radiusPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.radiusPt
	
	def clearRadius(self):
	
		for pt in self.radiusPt:
			pt.remove()
			
		self.radiusPt=[]
		
		self.fig.canvas.draw()
		
		return self.radiusPt
		
		
	def showFirstDataImg(self):
		
		self.embryo.updateFileList()
		
		fnImg=self.embryo.fnDatafolder
		if fnImg[-1]!='/':
			fnImg=fnImg+'/'
		fnImg=fnImg+self.embryo.fileList[0]
		
		img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
	
		self.showImg(img)
	
	def showImg(self,img):
		
		self.ax.imshow(img)
		self.fig.canvas.draw()
		
		return self.canvas
	
	def computeRadiusFromCoordinate(self,x,y):
		return np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
		
	def getMouseInput(self,event):
		
		#Check if click in axes
		if event.xdata==None:
			return
		
		#Left click to define center
		if event.button==1:
			
			self.center=[event.xdata,event.ydata]
			
			self.drawCenter()
			
			if len(self.radiusPt)>0:
				self.drawRadius()
				
		#Right click to define radius
		elif event.button==3:
			
			if len(self.centerPt)>0:
			
				self.radius=self.computeRadiusFromCoordinate(event.xdata,event.ydata)
				self.drawRadius()
			
		#Middle click to activate center marker	
		if event.button==2:
			
			self.drawCenterMarker()
				
		self.fig.canvas.draw()
		
		return
	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Generates matplotlib figure with (x,y) subplots

def makeSubplot(size,titles=None,tight=False,sup=None,proj=None,fig=None,show=True):
	
	#How many axes need to be created
	n_ax=size[0]*size[1]
	
	#Check if titles has right size
	if titles!=None:
		
		if len(titles)!=n_ax:
			print "Warning: len(titles)!=n_ax, will not print titles!"
			titles=None
	if proj!=None:
		if len(proj)!=n_ax:
			print "Warning: len(proj)!=n_ax, will not do projections!"
			proj=n_ax*[None]
	else:
		proj=n_ax*[None]
	
	#Creating figure
	if fig==None:
		fig=plt.figure()
		if show:
			fig.show()
	fig.set_tight_layout(tight)
	
	#Suptitle
	if sup!=None:
		fig.suptitle(sup)
	
	#Add axes
	axes=[]
	for i in range(n_ax):
		
		if proj[i]!=None:
			ax=fig.add_subplot(size[0],size[1],i+1,projection=proj[i])
		else:
			ax=fig.add_subplot(size[0],size[1],i+1)
		axes.append(ax)
		
		#Print titles 
		if titles!=None:
			ax.set_title(titles[i])
	
	#Draw
	plt.draw()
	
	#Return axes handle
	return fig,axes

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Adjust display range of imshow plots in list of axes

def adjustImshowRange(axes,vmin=None,vmax=None):
	
	#Loop through axes
	for ax in axes:
		
		#Grab plot
		implot=findImageArtist(ax,"Image")
		
		#Rescale
		if implot!=None:
			implot.set_clim(vmin,vmax)
		
	#Draw
	redraw(ax)
		
	return axes

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds artist with key in name and returns it

def findImageArtist(ax,key):
	
	c=ax.get_children()
	
	for x in c:
		if key in str(x):
			return x
			break
		
	return None	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Redraws axes's figure

def redraw(ax):
			
	ax.get_figure().canvas.draw()
	
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Plot timeseries all-in-one function

def plotTS(xvec,yvec,label='',title='',sup='',ax=None,color=None,linewidth=1,legend=True,linestyle='-'):
		
	if len(xvec)!=len(yvec):
		printWarning('len(xvec) != len (yvec). This could be due to incomplete simulation/analysis/pinning. Will not plot.')
		return None
	
	if ax==None:
		fig,axes = makeSubplot([1,1],titles=[title],sup=sup,tight=False)
		ax=axes[0]
	else:
		ax.set_title(title)
		
	ax.plot(xvec,yvec,color=color,label=label,linestyle=linestyle)
	
	if legend:
		ax.get_figure().subplots_adjust(right=0.7)
		ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	
	redraw(ax)
	
	return ax
		