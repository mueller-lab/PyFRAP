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
	
	"""PyFRAP simulation class. 
	
	Stores all important properties about how FRAP simulation is performed, such as:
	
		
	
	"""
	
	
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
		
		Initial condition modes are defined as:
		
			* 0: **ROI-based**: Mesh nodes get assigned the value of the first entry ``dataVec`` of 
			  the ROI covering them. Note: If a mesh node is covered by two ROIs, will assign the value
			  of the ROI that is last in embryo's ``ROIs`` list. See also 
			  :py:func:`pyfrp.modules.pyfrp_sim_module.applyROIBasedICs`.
			* 1: **Radial**: Concentrations get radially approximated around some center.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyRadialICs`.
			* 2: **Imperfect**: Interpolates ``ICimg`` onto mesh, fills nodes that are not covered
			  by image with ``embryo.analysis.concRim`` and then mimics imperfect bleaching via sigmoid
			  function in z-direction.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyImperfectICs`.
			* 3: **Inpterpolated**: Interpolates ``ICimg`` onto mesh, fills nodes that are not covered
			  by image with ``embryo.analysis.concRim``.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyInterpolatedICs`.
			* 4: **Ideal**: Ideal ICs, that is, single value inside and single value outside of bleached region.
			  See also :py:func:`pyfrp.modules.pyfrp_sim_module.applyIdealICs`, :py:func:`setValOut` 
			  and :py:func:`setBleachedROI`.
		
		.. note:: The default mode is **Interpolated** (``m=3``) and is highly recommended to obtain most realistic results.
		
		Args:
			m (int): Which mode to be used.
			
		Returns:
			int: Current initial condition mode used.
		
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
		
		"""Runs simulation.
		
		Checks if ROI indices are computed, if not, computes them. Then passes simulation
		object to :py:func:`pyfrp.modules.pyfrp_sim_module.simulateReactDiff`.
		
		Keyword Args:
			signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
			embCount (int): Counter of counter process if multiple datasets are simulated. 
			debug (bool): Print debugging messages and show debugging plots.
			showProgress (bool): Print out progress.
		
		Returns:
			pyfrp.subclasses.pyfrp_simulation.simulation: Updated simulation instance.
			
		
		"""
		
		if not self.embryo.checkROIIdxs()[1]:
			self.embryo.computeROIIdxs()
		
		pyfrp_sim_module.simulateReactDiff(self,signal=signal,embCount=embCount,showProgress=showProgress,debug=debug)
		return True
	
	def setMesh(self,m):
		
		"""Sets mesh to a new mesh object.
	
		Args:
			m (pyfrp.subclasses.pyfrp_mesh.mesh): PyFRAP mesh object.
			
		Returns:
			pyfrp.subclasses.pyfrp_mesh.mesh: Updated mesh instance.
			
		
		"""
		
		self.mesh=m
		return self.mesh
		
	def setICimg(self,img):
		
		"""Sets image for initial condition interpolation.
		
		Args:
			img (numpy.ndarray): A 2D image.
			
		Returns:
			numpy.ndarray: New ICimg.
		
		"""
		
		self.ICimg=img
		return self.ICimg
	
	def getICimg(self):
			
		"""Returns image for initial condition interpolation.
		
		Returns:
			numpy.ndarray: Current ICimg.
		
		"""
		
		
		return self.ICimg
	
	def showIC(self,ax=None,roi=None,nlevels=25,vmin=None,vmax=None):
		
		"""Plots initial conditions applied to mesh in 2D.
		
		If ``roi`` is given, will only plot initial conditions for nodes inside ROI, else 
		will plot initial condition for all nodes in mesh.
		
		.. note:: Simulation needs to be run first before this plotting function
		   can be used.
		
		.. image:: ../imgs/pyfrp_simulation/showIC.png
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI object.
			vmin (float): Overall minimum value to be displayed in plot.
			vmax (float): Overall maximum value to be displayed in plot.
			ax (matplotlib.axes): Axes used for plotting.
			nlevels (int): Number of contour levels to display.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		
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
		
		"""Plots image used for initial 2D.
		
		.. image:: ../imgs/pyfrp_simulation/showICimg.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes used for plotting.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		
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
			#plt_ICs=ax.imshow()
			
			#ax.autoscale(enable=True, axis='both', tight=True)
			plt.axis('equal')
			#cb=plt.colorbar(plt_ICs,orientation='horizontal',pad=0.05,shrink=0.9)
			
			plt.draw()
			
			return ax
			
		else:
			printWarning("ICimg is not analyzed yet. Run data analysis first.")
			return None
	
	def computeInterpolatedICImg(self):
		
		"""Computes interpolation of initial condition image.
		
		Interpolation is done as in :py:func:`pyfrp.modules.pyfrp_sim_module.applyInterpolatedICs`.
		
		Returns:
			tuple: Tuple containing:
			
				* xInt (numpy.ndarray): Meshgrid x-coordinates.
				* yInt (numpy.ndarray): Meshgrid y-coordinates.
				* f (numpy.ndarray): Interpolated image.
	
		"""
		
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
		
		"""Interpolates ICs back onto 2D image.
		
		Uses ``scipy.interpolate.griddata``, see also http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.griddata.html
		
		If ``roi`` is specified, will only interpolate nodes of this ROI. 
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			
		Returns:
			tuple: Tuple containing:
			
				* X (numpy.ndarray): Meshgrid x-coordinates.
				* Y (numpy.ndarray): Meshgrid y-coordinates.
				* interpIC (numpy.ndarray): Interpolated ICs.
		
		"""
		
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
		
		"""Shows ICs interpolated back onto 2D image.
		
		If ``roi`` is specified, will only interpolate nodes of this ROI. 
		
		See also :py:func:`computeInterpolatedIC`.
		
		.. image:: ../imgs/pyfrp_simulation/showInterpolatedIC.png
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			ax (matplotlib.axes): Axes to be used for plotting.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		X,Y,interpIC=self.computeInterpolatedIC(roi=roi)
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Interpolated IC"],sup="simulation")
			ax=axes[0]
			
		ax.imshow(interpIC)
		
		return ax
	
	def showInterpolatedICImg(self,ax=None):
		
		"""Shows interpolation of initial condition image.
		
		See also :py:func:`computeInterpolatedICImg`.
		
		.. image:: ../imgs/pyfrp_simulation/showInterpolatedICimg.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be used for plotting.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
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
		
		"""Shows initial image, its interpolation, the resulting initial 
		condition and its interpolation back onto an image.
		
		See also :py:func:`showICimg`, :py:func:`showInterpolatedICImg`, 
		:py:func:`showIC`, :py:func:`showInterpolatedIC`.
		
		Will create new axes if necessary.
		
		.. warning:: Some images might be flipped due to plotting functions. Will be fixed in future version.
		
		.. image:: ../imgs/pyfrp_simulation/ICcompare.png
		
		Keyword Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			axes (matplotlib.axes): List of axes of length 4.
			
		Returns:
			list: List of axes.
		
		"""
	
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
		
		r"""Generates time vector that is optimal to fit 
		experiments with expected diffusion coefficients up
		to ``maxDExpectedPx``.
		
		Basically computes how long a simulation needs to run in
		seconds to capture the dynamics of an experiment with diffusion
		coefficient of ``maxDExpectedPx``. Does this by setting end time point to
		
		.. math:: t_{\mathrm{end,sim}} = \frac{D_{\mathrm{max. exp.}}}{D_{\mathrm{sim}}} t_{\mathrm{end,data}}
		
		.. note:: Keeps time scaling.
		
		Args:
			maxDExpectedPx (float): Maximum expected diffusion coefficient.
			
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		wasLog=self.isLogTimeScale()
			
		self.tvecSim=np.linspace(self.embryo.tStart,maxDExpectedPx/self.D*self.embryo.tEnd,self.stepsSim)
		
		if wasLog:
			self.toLogTimeScale()
		
		return self.tvecSim
		
	def toLogTimeScale(self,spacer=1E-10):
		
		"""Converts time vector for simulation to logarithmic scale.
		
		Keyword Args:
			spacer (float): Small offset to avoid log(0).
			
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		self.tvecSim=self.tvecSim[0]+np.logspace(np.log10(spacer+self.tvecSim[0]), np.log10(self.tvecSim[-1]), self.stepsSim)-spacer
		return self.tvecSim
		
	def toLinearTimeScale(self):
		
		"""Converts time vector for simulation to linear scale.
	
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		self.tvecSim=np.linspace(self.tvecSim[0],self.tvecSim[-1],self.stepsSim)	
		return self.tvecSim
	
	def toDefaultTvec(self):
		
		"""Sets time vector for simulation to default range.
		
		Default range is given by ``tStart`` and ``tEnd`` in ``embryo`` object
		and is linearly scaled.
		
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		self.tvecSim=np.linspace(self.embryo.tStart,self.embryo.tEnd,self.stepsSim)
		return self.tvecSim
	
	def updateTvec(self):
		
		"""Updates time vector for simulation to match 
		experiment start and end time.
		
		Does not change scaling of time vector.
		
		Returns:
			numpy.ndarray: New simulation time vector.
		
		"""
		
		if (self.tvecSim[1]-self.tvecSim[0])==(self.tvecSim[-1]-self.tvecSim[-2]):
			self.toLinearTimeScale()
		else:
			self.toLogTimeScale()
		return self.tvecSim
	
	def setTimesteps(self,n):
		
		"""Sets number of simulation time steps and updates
		time vector.
		
		Args:
			n (int): New number of time steps.
			
		Returns:
			int: New number of time steps.
			
		"""
		
		self.stepsSim=int(n)
		self.updateTvec()
		return self.stepsSim
	
	def setD(self,D):
		
		"""Sets diffusion coefficient used for simulation.
		
		Args:
			D (float): New diffusion coefficient in :math:`\mu\mathrm{m}^2/s`.
			
		Returns:
			float: New diffusion coefficient in :math:`\mathrm{px}^2/s`.
			
		"""
		
		self.D=D
		return self.D
	
	def getD(self):
		
		"""Returns current diffusion coefficient used for simulation.
		
		Returns:
			float: Current diffusion coefficient in :math:`\mathrm{px}^2/s`.
			
		"""
		
		return self.D
	
	def setProd(self,prod):
		
		"""Sets production rate used for simulation.
		
		Args:
			prod (float): New production rate in :math:`1/s`.
			
		Returns:
			float: New production rate in :math:`1/s`.
			
		"""
		
		self.prod=prod
		return self.prod
	
	def getProd(self):
		
		"""Returns production rate used for simulation.
		
		Returns:
			float: Current production rate in :math:`1/s`.
			
		"""
		
		return self.prod
	
	def setDegr(self,degr):
		
		"""Sets degradation rate used for simulation.
		
		Args:
			prod (float): New degradation rate in :math:`1/[c]s`.
			
		Returns:
			float: New degradation rate in :math:`1/[c]s`.
			
		"""
		
		self.degr=degr
		return self.degr
	
	def getDegr(self):
		
		"""Returns degradation rate used for simulation.
			
		Returns:
			float: Current degradation rate in :math:`1/[c]s`.
			
		"""
		
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
		
		.. image:: ../imgs/pyfrp_simulation/ICstack.png
		
		
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
		
		
		
		
	