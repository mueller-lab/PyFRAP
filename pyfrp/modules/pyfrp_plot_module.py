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

"""Plotting module for PyFRAP toolbox. 

Contains functions and classes that are often used by PyFRAP toolbox and simplify plot creation and management.

"""

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

class FRAPBoundarySelector():
	
	"""Simple GUI widget to select circular FRAP boundaries.
	
	Has useful center marker that is activatable that helps finding the 
	center of the image.
	
	Mouse Input:
		
		* Left: Set center.
		* Right: Set Radius.
		* Middle: Activate center marker.
		
	Keyboard Input:
	
		* Left Arrow: Move center to the left.
		* Right Arrow: Move center to the right.
		* Up Arrow: Move center upwards.
		* Down Arrow: Move center downwards.
		* Control + Up Arrow: Increase circle radius.
		* Control + Down Arrow: Decrease circle radius.
	
	.. note:: If ``embryo`` is given at initiation, will use first image specified
	   in embryo's ``fileList`` as background image.
	   
	Example Usage:
	
	>>> sel=FRAPBoundarySelector(fn="path/to/img/file")

	Use mouse/keyboard to define circular boundary.
	
	>>> center,radius=sel.getResults()
	
	Keyword Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): PyFRAP embryo instance.
		fn (str): Filepath to image file taken for boundary selection.
		
	"""
	
	def __init__(self,embryo=None,fn=None):
		
		#Passing embryo to class
		self.embryo=embryo	
		self.fn=fn
		
		#Img
		self.img=None
		
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
		
		#Close everthing down if faulty input
		self.checkInput()
			
		#Plot image if existent
		self.showFirstDataImg()
		plt.show()
	
	def checkInput(self):
		
		"""Checks if at least one of the two keyword arguments, ``embryo`` or ``fn``, is given.
		
		If not, prints error message and closes down widget.
		"""
		
		if self.fn==None and self.embryo==None:
			printError("No Embryo or fn defined. Going to exit.")
			plt.close(self.fig)
			return
	
	def createCanvas(self):
		
		"""Creates figure and canvas used for plotting."""
			
		h=500/self.dpi
		v=500/self.dpi
		
		self.fig = plt.figure()
		self.fig.show()
		self.fig.set_size_inches(h,v,forward=True)
		
		self.canvas=self.fig.canvas
		
		self.ax = self.fig.add_subplot(111)
		
		if self.embryo!=None:
			self.ax.set_xlim([0, self.embryo.dataResPx])
			self.ax.set_ylim([0, self.embryo.dataResPx])
			
		#Connect to mouseclick event
		self.fig.canvas.mpl_connect('close_event', self.closeFigure)
		self.canvas.mpl_connect('button_press_event', self.getMouseInput)
		self.canvas.mpl_connect('key_press_event', self.keyPressed)
		
		self.canvas.draw()
		
		return 
	
	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		elif event.key=='ctrl+up':
			self.increaseRadius()
		elif event.key=='ctrl+down':
			self.decreaseRadius()
			
	def moveLeft(self):
		
		"""Moves center 1 px to the left."""
		
		if self.center!=None:
			self.center=[self.center[0]-1,self.center[1]]
			self.drawCenter()
		
		
	def moveRight(self):
		
		"""Moves center 1 px to the right."""
		
		if self.center!=None:
			self.center=[self.center[0]+1,self.center[1]]
			self.redraw()
	
	def moveUp(self):
		
		"""Moves center 1 px up."""
		
		if self.center!=None:
			self.center=[self.center[0],self.center[1]+1]
			self.redraw()
	
	def moveDown(self):
		
		"""Moves center 1 px down."""
		
		if self.center!=None:
			self.center=[self.center[0],self.center[1]-1]
			self.redraw()
	
	def redraw(self):
		
		"""Redraws both center and radius if available."""
		
		if self.center!=None:
			self.drawCenter()
		if self.radius!=None:
			self.drawRadius()
			
	def increaseRadius(self):
		
		"""Increases radius by 1 px."""
		
		if self.radius!=None:
			self.radius=self.radius+1
			self.redraw()
			
	def decreaseRadius(self):
		
		"""Decreases radius by 1 px."""
		
		if self.radius!=None:
			self.radius=self.radius-1
			self.redraw()
			
	def closeFigure(self,event):
		
		"""Returns center and radius at close ``event``."""
		
		return self.center,self.radius
		
	def getEmbryo(self):	
		
		"""Returns ``embryo`` object if given.``"""
		
		return self.embryo
	
	def getResults(self):
		
		"""Returns center and radius."""
		
		return self.center,self.radius
	
	def drawCenterMarker(self):
		
		"""Draws a yellow marker in center of the image, making it 
		easier to find image center when selecting center of boundary."""
		
		centerImg=[self.img.shape[0]/2.,self.img.shape[1]/2.]
		
		if len(self.centerMarker)>0:
			self.clearCenterMarker()
		else:	
			pt=ptc.Circle(centerImg,radius=3,fill=True,color='y')
			self.centerMarker.append(self.ax.add_patch(pt))
			
		self.fig.canvas.draw()
		
		return self.centerMarker
		
	def clearCenterMarker(self):
		
		"""Removes center maker from canvas."""
	
		for pt in self.centerMarker:
			pt.remove()
			
		self.fig.canvas.draw()	
			
		self.centerMarker=[]
	
	def drawCenter(self):
		
		"""Draws a red marker at selected center on canvas."""
		
		if len(self.centerPt)>0:
			self.clearCenter()
	
		pt=ptc.Circle(self.center,radius=3,fill=True,color='r')
		self.centerPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.centerPt
	
	def clearCenter(self):
	
		"""Removes center marker from canvas."""
		
		for pt in self.centerPt:
			pt.remove()
			
		self.centerPt=[]
		
		self.fig.canvas.draw()
		
	def drawRadius(self):
		
		"""Draws a red circle around selected center with selected radius on canvas."""
		
		if len(self.radiusPt)>0:
			self.clearRadius()
			
		pt=ptc.Circle(self.center,radius=self.radius,fill=False,color='r')
		self.radiusPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.radiusPt
	
	def clearRadius(self):
		
		"""Removes circle from canvas."""
		
		for pt in self.radiusPt:
			pt.remove()
			
		self.radiusPt=[]
		
		self.fig.canvas.draw()
		
		return self.radiusPt
		
		
	def showFirstDataImg(self):
		
		"""Shows either first data image defined in ``embryo.fileList`` or 
		image specified by ``fn``.
		
		.. note:: If both are given, will use the embryo option.
		
		"""
		
		if self.embryo!=None:
			self.embryo.updateFileList()
		
			fnImg=self.embryo.fnDatafolder
			if fnImg[-1]!='/':
				fnImg=fnImg+'/'
			fnImg=fnImg+self.embryo.fileList[0]
			
			self.img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
		
		elif self.fn!=None:
			self.img=pyfrp_img_module.loadImg(self.fn,'uint16')
			self.ax.set_xlim([0, self.img.shape[0]])
			self.ax.set_ylim([0, self.img.shape[1]])
			
		self.showImg(self.img)
	
	def showImg(self,img):
		
		"""Shows image on canvas.
		
		Args:
			img (numpy.ndarray): Image to be shown.
		
		"""
		
		self.ax.imshow(img)
		self.fig.canvas.draw()
		
		return self.canvas
	
	def computeRadiusFromCoordinate(self,x,y):
		
		"""Computes radius from given cordinate ``(x,y)``.
		"""
		
		return np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
		
	def getMouseInput(self,event):
		
		"""Directs mouse input to the right actions."""
		
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
	
def makeSubplot(size,titles=None,tight=False,sup=None,proj=None,fig=None,show=True):
	
	"""Generates matplotlib figure with (x,y) subplots.
	
	.. note:: List of ``titles`` needs to be same size as desired number of axes.
	   Otherwise will turn off titles.
	
	.. note:: List of ``proj`` needs to be same size as desired number of axes.
	   Otherwise will turn off projections.
	   
	Example:
	
	>>> makeSubplot([2,2],titles=["Axes 1", "Axes 2", "Axes 3", "Axes 4"],proj=[None,None,'3d',None])
	
	will result in
	
	.. image:: ../imgs/pyfrp_plot_module/makeSubplot.png
	
	Args:
		size (list): Size of subplot arrangement.
	
	Keyword Args:
		titles (list): List of axes titles.
		tight (bool): Use tight layout.
		sup (str): Figure title.
		proj (list): List of projections.
		fig (matplotlib.figure): Figure used for axes.
		show (bool): Show figure right away.
	
	Returns:
		tuple: Tuple containing:
			
			* fig (matplotlib.figure): Figure.
			* axes (list): List of Matplotlib axes.
	
	"""
	
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

def adjustImshowRange(axes,vmin=None,vmax=None):
	
	"""Adjust display range of ``matplotlib.pyplot.imshow`` plots in 
	list of axes.
	
	Finds first image artist in each axes in ``axes`` list and then
	sets display range to ``[vmin,vmax]``.
	
	Args:
		axes (list): List of matplotlib axes.
		
	Keyword Args:
		vmin (float): Minimum value of display range.
		vmax (float): Maximum value of display range.
		
	Returns:
		 list: Updated list of matplotlib axes.
	
	"""
	
	#Loop through axes
	for ax in axes:
		
		#Grab plot
		implot=findArtist(ax,"Image")
		
		#Rescale
		if implot!=None:
			implot.set_clim(vmin,vmax)
		
	#Draw
	redraw(ax)
		
	return axes

def findArtist(ax,key):
	
	"""Finds ``matplotlib.artist`` which name contains ``key``.
	
	.. note:: Will stop searching after first artist is found.
	
	Will return ``None`` if no artist can be found.
	
	Args:
		ax (matplotlib.axes): Matplotlib axes.
		key (str): Key used for search.
		
	Returns:
		matplotlib.artist: Matplotlib artist.
	
	"""
	
	c=ax.get_children()
	
	for x in c:
		if key in str(x):
			return x
			break
		
	return None	

def redraw(ax):
	
	"""Redraws axes's figure's canvas.
	
	Makes sure that current axes content is visible.
	
	Args:
		ax (matplotlib.axes): Matplotlib axes.
		
	Returns:
		matplotlib.axes: Matplotlib axes 
	
	"""
			
	ax.get_figure().canvas.draw()
	
	return ax

def plotTS(xvec,yvec,label='',title='',sup='',ax=None,color=None,linewidth=1,legend=True,linestyle='-'):
	
	"""Plot timeseries all-in-one function.
	
	Args:
		xvec (numpy.ndarray): x-data to be plotted.
		yvec (numpy.ndarray): y-data to be plotted.
	
	Keyword Args:
		ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
		color (str): Color of plot.
		linestyle (str): Linestyle of plot.
		linewidth (float): Linewidth of plot.
		legend (bool): Show legend.
		sup (str): Figure title.
		title (str): Axes title.
		label (str): Label for legend.
	
	Returns:
		matplotlib.axes: Axes used for plotting.
	
	"""
	
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

	