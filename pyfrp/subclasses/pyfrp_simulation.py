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

#Simulation class for PyFRAP toolbox, including following classes:

#(1) simulation

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np
import scipy.interpolate as interp 

#PyFRAP classes
import pyfrp_mesh

#PyFRAP modules
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_sim_module
from pyfrp.modules.pyfrp_term_module import *

#Plotting
import matplotlib.pyplot as plt

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================

class simulation(object):

	def __init__(self,embryo):
		
		#Naming/ID
		self.embryo=embryo
		
		#IC image
		self.ICimg=None
		
		self.restoreDefaults()
		
	def restoreDefaults(self):
		
		#PDE specific
		self.D=50.
		self.prod=0.0
		self.degr=0.0
		
		#Mesh specific
		self.mesh=self.setMesh(pyfrp_mesh.mesh(self))
		
		#IC specific
		self.ICmode=3
		self.IC=None
		
		#Only necessary if ICmode=4 (ideal ICs)
		self.bleachedROI=None
		self.valOut=None
		
		#Time specific
		self.stepsSim=3000
		self.tvecSim=np.linspace(self.embryo.tStart,self.embryo.tEnd,self.stepsSim)
		
		#Integration specific (deprecated)
		#self.avgMode=0
		#self.addRimSim=0
		#self.intSteps=512
		#self.integrationMethod=2
		
		#Debugging (deprecated)
		#self.debug=0
		
	def setICMode(self,m):
		
		"""Sets the mode of initial conditions.
		"""
		
		if m not in range(5):
			printError("ICmode = " +  m + " is not defined. Not going to change ICmode" )
			return self.ICmode
			
		self.ICmode=m
		
		if m==4:
			if self.bleachedROI==None:
				printWarning("bleachedROI is not set yet. This might lead to problems later.")
			if self.valOut==None:
				printWarning("valOut is not set yet. This might lead to problems later.")
		
		return self.ICmode
	
	def setBleachedROI(self,r):
		
		"""Sets bleached ROI that is used when ideal ICs (ICmode=4) is selected.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to be set bleached ROI.
			
		"""
		
		self.bleachedROI=r
		
		return self.bleachedROI
	
	def setValOut(self,v):
		
		"""Sets valOut that is used when ideal ICs (ICmode=4) is selected.
		
		Args:
			v (float): Value that is to assigned outside of bleachedROI.
		
		"""
		
		self.valOut=v
		
		return self.valOut
	
	def run(self,signal=None,embCount=None,showProgress=True,debug=False):
		
		if not self.embryo.checkROIIdxs()[1]:
			self.embryo.computeROIIdxs()
		
		pyfrp_sim_module.simulateReactDiff(self,signal=signal,embCount=embCount,showProgress=showProgress,debug=debug)
		return True
	
	def setMesh(self,m):
		self.mesh=m
		return self.mesh
		
	def setICimg(self,img):
		self.ICimg=img
		return self.ICimg
	
	def getICimg(self):
		return self.ICimg
	
	def showIC(self,ax=None,roi=None,nlevels=25,vmin=None,vmax=None):
		
		if self.IC!=None:
			
			if ax==None:
				fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["IC"],sup="",tight=False)
				ax=axes[0]
			
			if roi==None:
		
				x,y,z=self.mesh.mesh.getCellCenters()
			
				if vmin==None:
					vmin=min(self.IC)
				if vmax==None:
					vmax=max(self.IC)
			
				levels=np.linspace(vmin,1.01*vmax,nlevels)
		
				ax.tricontourf(x,y,self.IC,vmin=vmin,vmax=vmax,levels=levels)
				ax.autoscale(enable=True, axis='both', tight=True)
			
				ax.get_figure().canvas.draw()
		
			else:
				ax=roi.plotSolutionVariable(self.IC,ax=ax,nlevels=nlevels,vmin=vmin,vmax=vmax)
			
			return ax
			
		else:
			printWarning("IC is not generated yet. Run simulation first.")
			return None	
		
	def showICimg(self,ax=None):
		
		if self.ICimg!=None:
			
			if ax==None:
				fig,axes = pyfrp_plot_module.makeSubplot([1,1],sup="",tight=False)
				ax=axes[0]
				
			res=self.ICimg.shape[0]
			if 'quad' in self.embryo.analysis.process.keys():
				X,Y=np.meshgrid(np.arange(res,2*res),np.arange(res,2*res))
			else:
				X,Y=np.meshgrid(np.arange(res),np.arange(res))
			
			plt_ICs=ax.contourf(X,Y,self.ICimg)
			
			ax.autoscale(enable=True, axis='both', tight=True)
			plt.axis('equal')
			#cb=plt.colorbar(plt_ICs,orientation='horizontal',pad=0.05,shrink=0.9)
			
			plt.draw()
			
			return ax
			
		else:
			printWarning("ICimg is not analyzed yet. Run data analysis first.")
			return None
	
	def computeInterpolatedICImg(self):
		
		#Get image resolution and center of geometry
		res=self.ICimg.shape[0]
		center=self.embryo.geometry.getCenter()
		
		#Define x/y coordinates of interpolation
		if 'quad' in self.embryo.analysis.process.keys():
			#Shift everything by center to fit with the mesh
			xInt = np.arange(center[0]+1, center[0]+res+1, 1)
			yInt = np.arange(center[1]+1, center[1]+res+1, 1)		
		else:
			xInt = np.arange(1, res+1, 1)
			yInt = np.arange(1, res+1, 1)
		
		#Generate interpolation function
		f=interp.RectBivariateSpline(xInt, yInt, self.ICimg, bbox=[None, None, None, None], kx=3, ky=3, s=0)
		
		return xInt, yInt, f
	
	def computeInterpolatedIC(self,roi=None):
		
		#Get image resolution and center of geometry
		res=self.ICimg.shape[0]
		
		#Build Empty Img
		X,Y=np.meshgrid(np.arange(res),np.arange(res))
		
		#Get cellcenters
		x,y,z=self.mesh.mesh.getCellCenters()
		
		if roi!=None:
			xInt=x[roi.meshIdx]
			yInt=y[roi.meshIdx]
			val=self.IC[roi.meshIdx]
		else:
			xInt=x
			yInt=y
			val=self.IC
		
		interpIC=interp.griddata((xInt,yInt),val,(X,Y),method='linear')
		
		return X,Y,interpIC
		
	def showInterpolatedIC(self,ax=None,roi=None):
		
		X,Y,interpIC=self.computeInterpolatedIC(roi=roi)
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Interpolated IC"],sup="simulation")
			ax=axes[0]
			
		ax.imshow(interpIC)
		
		return ax
	
	def showInterpolatedICImg(self,ax=None):
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Interpolated Image"],sup="simulation")
			ax=axes[0]
			
		xInt, yInt, f=self.computeInterpolatedICImg()	
		
		#print np.shape(xInt)
		
		#raw_input()
		X,Y=np.meshgrid(xInt,yInt)
		imgInt=np.zeros(np.shape(X))
		
		for i in range(np.shape(X)[0]):
			for j in range(np.shape(Y)[0]):
				imgInt[i,j]=f(X[i,j],Y[i,j])
				
		ax.imshow(imgInt)	
		
		return ax
	
	def compareICInterpolation(self,axes=None,roi=None):
	
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([2,2],titles=["Original Image","Interpolated Image","IC","Reinterpolated IC"],sup="simulation")
		
		self.showICimg(ax=axes[0])
		self.showInterpolatedICImg(ax=axes[1])
		self.showIC(ax=axes[2],roi=roi)
		self.showInterpolatedIC(ax=axes[3],roi=roi)
		
		for ax in axes:
			pyfrp_plot_module.redraw(ax)
		
		return axes	
		
	def getOptTvecSim(self,maxDExpectedPx):
		self.tvecSim=np.linspace(self.embryo.tStart,maxDExpectedPx/self.D*self.embryo.tEnd,self.stepsSim)
		return self.tvecSim
		
	def toLogTimeScale(self,spacer=1E-10):
		self.tvecSim=self.tvecSim[0]+np.logspace(np.log10(spacer+self.tvecSim[0]), np.log10(self.tvecSim[-1]), self.stepsSim)-spacer
		return self.tvecSim
		
	def toLinearTimeScale(self):
		self.tvecSim=np.linspace(self.tvecSim[0],self.tvecSim[-1],self.stepsSim)	
		return self.tvecSim
	
	def toDefaultTvec(self):
		self.tvecSim=np.linspace(self.embryo.tStart,self.embryo.tEnd,self.stepsSim)
		return self.tvecSim
	
	def updateTvec(self):
		if (self.tvecSim[1]-self.tvecSim[0])==(self.tvecSim[-1]-self.tvecSim[-2]):
			self.toLinearTimeScale()
		else:
			self.toLogTimeScale()
		return self.tvecSim
	
	def setTimesteps(self,n):
		self.stepsSim=int(n)
		self.updateTvec()
		return self.stepsSim
	
	def setD(self,D):
		self.D=D
		return self.D
	
	def getD(self):
		return self.D
	
	def setProd(self,prod):
		self.prod=prod
		return self.prod
	
	def getProd(self):
		return self.prod
	
	def setDegr(self,degr):
		self.degr=degr
		return self.degr
	
	def getDegr(self):
		return self.degr
	
	def isLogTimeScale(self):
		
		"""Returns if time spacing of simulation is logarithmic.
		
		Returns:
			bool: Time spacing is logarithmic.
		
		"""
		
		#Note: We round here to 5 decimals since diff seems to make numerical mistakes.
		
		return not round(np.diff(self.tvecSim)[0],5)==round(np.diff(self.tvecSim)[-1],5)
	
	def plotICStack(self,ROIs,withGeometry=True,vmin=None,vmax=None,ax=None,colorbar=False):
		
		"""Plots a stack of the initial conditions in a given list of ROIs.
		
		Will automatically compute the direction in which ROI lies in the 3D space and
		reduce the ROI into this plane for contour plot.
		
		If ``vmin=None`` or ``vmax=None``, will compute overall maximum and minimum values
		over all ROIs.
		
		Args:
			phi (fipy.CellVariable): Simulation solution variable (or numpy array).
			ROIs (list): List of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` objects.
			
		Keyword Args:
			withGeometry (bool): Show geometry inside plot.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			ax (matplotlib.axes): Axes used for plotting.
			colorbar (bool): Display color bar.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		ax = self.plotSolStack(self.IC,ROIs,withGeometry=withGeometry,vmin=vmin,vmax=vmax,ax=ax,colorbar=colorbar)
		
		return ax 
	
	
		
	def plotSolStack(self,phi,ROIs,withGeometry=True,vmin=None,vmax=None,ax=None,colorbar=False):
			
		"""Plots a stack of the solution variable in a given list of ROIs.
		
		Will automatically compute the direction in which ROI lies in the 3D space and
		reduce the ROI into this plane for contour plot.
		
		If ``vmin=None`` or ``vmax=None``, will compute overall maximum and minimum values
		over all ROIs.
		
		Args:
			phi (fipy.CellVariable): Simulation solution variable (or numpy array).
			ROIs (list): List of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` objects.
			
		Keyword Args:
			withGeometry (bool): Show geometry inside plot.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			ax (matplotlib.axes): Axes used for plotting.
			colorbar (bool): Display color bar.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
			
			
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Simulation IC stack"],proj=['3d'])
			ax=axes[0]
		
		#Plot geometry
		if withGeometry:
			#self.embryo.geometry.updateGeoFile()
			ax=self.embryo.geometry.plotGeometry(ax=ax)
		
		#Find vmin/vmax over all ROIs
		vminNew=[]
		vmaxNew=[]
		for r in ROIs:
			vminNew.append(min(phi[r.meshIdx]))
			vmaxNew.append(max(phi[r.meshIdx]))
		
		if vmin==None:
			vmin=min(vminNew)	
		if vmax==None:
			vmax=min(vmaxNew)
			
		for r in ROIs:
			plane=r.getMaxExtendPlane()
			zs=r.getPlaneMidCoordinate()
			zdir=r.getOrthogonal2Plane()
			
			ax=r.plotSolutionVariable(phi,ax=ax,vmin=vmin,vmax=vmax,plane=plane,zs=zs,zdir=zdir,colorbar=colorbar)
				
		return ax
		
		
		
		
	