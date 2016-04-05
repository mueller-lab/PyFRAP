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

#Simulaton module for PyFRAP toolbox, including following functions:

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#PDE Toolbox
from fipy import *

#Numpy/Scipy
import numpy as np
import scipy.interpolate as interp 
import scipy.ndimage.interpolation as ndi

#matplotlib
import matplotlib.pyplot as plt

#Misc
import time
import sys

#PyFRAP Modules
import pyfrp_plot_module 
import pyfrp_integration_module
import pyfrp_misc_module
from pyfrp_term_module import *
import pyfrp_idx_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def simulateReactDiff(simulation,signal=None,embCount=None,showProgress=True,debug=False):
	
	r"""Simulates reaction diffusion equation goverining FRAP experiment.
	
	Performs the following steps:
		
		* Resets ``simVecs`` of all ROIs.
		* If not generated yet, generates mesh (should never be the case!)
		* Initializes PDE with Neumann boundary conditions resulting in the problem:
		
		  .. math::
		     \partial_t c = D \nabla^2 c - k_1 c + k_2,
	
		  where :math:`k_1` is the degradation rate and :math:`k_2` the production rate.
		
		* Applies initial conditions defined in ``simulation.ICmode``.
		* Simulates FRAP experimment.
		
	Args: 
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Keyword Args:
		signal (PyQt4.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
		embCount (int): Counter of counter process if multiple datasets are analyzed. 
		debug (bool): Print final debugging messages and show debugging plots.
		showProgress (bool): Show simulation progress. 
		
	Returns: 
		pyfrp.subclasses.pyfrp_simulation.simulation: Updated simulation object.
	"""
	

	#Stepping and timescale
	timeStepDuration = simulation.tvecSim[1]-simulation.tvecSim[0]
	
	#Reset simulation vecs
	for r in simulation.embryo.ROIs:
		r.resetSimVec()
	
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	print "Starting simulation"
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	startTimeTotal=time.clock()
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Mesh Generation
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	startTimeMesh=time.clock()
	
	if simulation.mesh.mesh==None:
		printWarning('No mesh has been generated yet!')
		a=raw_input('Do you want to generate a mesh now?[Y/N]')
		if a=='Y':
			simulation.mesh.genMesh()
			print "Mesh created in", time.clock()-startTimeMesh
		else:
			print 'Cannot run simulation without mesh, will abort.'
			return simulation
		
		timeMesh=time.clock()-startTimeMesh
		
	print "Mesh created after", time.clock()-startTimeTotal

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Initialization of PDE
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Create solution variable
	
	phi = CellVariable(name = "solution variable",mesh = simulation.mesh.mesh,value = 0.) 
	

	#Apply initial conditions
	if simulation.ICmode=='ideal':
		for r in simulation.embryo.ROIs:
			phi.value[r.meshIdx]=r.dataVec[0]
		
	elif simulation.ICmode=='radial':
		phi = applyRadialICs(phi,simulation,debug=debug)
	
	elif simulation.ICmode==2:
		phi = CellVariable(name = "solution variable",mesh = embryo.mesh,value = 0.) 
		phi=mimic_imperfect_bleaching(embryo.slice_height_px[0],embryo.cylinder_height_px,embryo.im_reg_ICs,embryo.conc_rim,embryo.rim,embryo.add_rim_from_radius,embryo.mesh,phi,embryo.side_length_bleached_px,embryo.center_embr_px,1)
		
	elif simulation.ICmode==3:
		phi=applyInterpolatedICs(phi,simulation,debug=False)
		
	#Remember ICs
	simulation.IC=np.asarray(phi.value).copy()
	

	#Defining Type of equation
	eq = TransientTerm() == DiffusionTerm(coeff=simulation.D)+simulation.prod-simulation.degr*phi

	#Defining BCs
	#Note: BCs are Neumann boundaries by default 
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Calculating initial concentrations 
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	for r in simulation.embryo.ROIs:
		r.getSimConc(phi,append=True)
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Solving PDE
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Keeping track of time
	startTimeSim=time.clock()
	
	avgTime=0
	stepTime=0

	#mySolver = LinearLUSolver(iterations=10, tolerance=5e-1)
	mySolver = LinearPCGSolver(tolerance=1e-15,iterations=5000)
	
	#solvers.DefaultSolver
	for step in range(simulation.stepsSim-1):
		
		#Compute timestep duration 
		timeStepDuration=simulation.tvecSim[step+1]-simulation.tvecSim[step]
		
		#Solve PDE in this Step
		stepStart=time.clock()
		eq.solve(var=phi,dt=timeStepDuration,solver=mySolver)
		stepTime=stepTime+(time.clock()-stepStart)
				
		#Compute concentration
		avgStart=time.clock()
		
		for r in simulation.embryo.ROIs:
			r.getSimConc(phi,append=True)
		
		avgTime=avgTime+(time.clock()-avgStart)
			
		#Print Progress
		if showProgress:
			currPerc=int(100*step/float(simulation.stepsSim))
			
			if signal==None:
				sys.stdout.write("\r%d%%" %currPerc)  
				sys.stdout.flush()
			else:	
				if embCount==None:
					signal.emit(currPerc)
				else:
					signal.emit(currPerc,embCount)
			

	print "Step time: ", stepTime, " in %:", stepTime/(time.clock()-startTimeSim)*100
	print "Avg time: ", avgTime, " in %:", avgTime/(time.clock()-startTimeSim)*100
	print "Simulation done after", time.clock()-startTimeTotal
	
	return simulation

def applyRadialICs(phi,simulation,radSteps=15,debug=False):
	
	"""Applies radially averaged image data to solution variable as IC.
	
	.. note:: Will use ``embryo.geometry.center`` as center circle and the maximum 
	   distant pixel from the center as maximum radius.
	
	Args: 
		phi (fipy.CellVariable): PDE solution variable.
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Keyword Args:
		radSteps (int): Number of radial levels.
		debug (bool): Print debugging messages.
		
	
	
	"""
	
	#Adjust center so histogram works for 'quad'
	if 'quad' in simulation.embryo.analysis.process.keys():
		center=[0,0]
	else:
		center=simulation.ICimg,simulation.embryo.geometry.getCenter()
	
	#Compute radial histogram of IC image
	maxR=pyfrp_img_module.dist(center,[simulation.ICimg.shape[0],simulation.ICimg.shape[0]])
	bins,binsMid,histY,binY=pyfrp_img_module.radialImgHist(simulation.ICimg,nbins=radSteps,byMean=True,maxR=maxR)
	
	#Set center to actual center of geometry
	center=simulation.ICimg,simulation.embryo.geometry.getCenter()
	
	#Apply value of most outer bin to all nodes
	phi.setValue(binY[-1])
	
	#Loop through all bins and apply values from outside to inside
	binsY=binsY.reverse()
	bins=bins.reverse()
	for i in range(len(binY)):
		
		phi.setValue(binY[i], where=(x-center[0])**2+(y-center[1])**2 < bins[i]**2)
			
		if debug:
			print "Applied concentration", binY[i], " to all nodes with radius <", bins[i] 
	
	return phi
		

def applyInterpolatedICs(phi,simulation,matchWithMaster=True,debug=False):
	
	"""Interpolates initial conditions onto mesh.
	
	Uses a bivarariate spline interpolation (http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.interpolate.RectBivariateSpline.html)
	to generate an interpolation function of the IC image. Then applies interpolated values to solution variable ``phi`` if mesh nodes are inside 
	image and masterROI. If not, will apply rim concentration.
	
	.. note:: If no rim concentration has been calculated (for example through running the data analysis)
	   applyInterpolatedICs will try to compute ``concRim`` by itself. For this it will take the mean concentration outside of bleached square but inside
	   ``masterROI``. 
	   
	.. note:: The bleached square used here is not defined as a ``ROI`` object here, but rather through the properties 
	   ``embryo.sideLengthBleachedPx`` and ``embryo.offsetBleachedPx``. This might change in future versions.
	
	Args: 
		phi (fipy.CellVariable): PDE solution variable.
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	Keyword Args:
		matchWithMaster (bool): Match interpolation indices with ``masterROI`` indices.
		debug (bool): Print debugging messages.
	
	Returns:
		fipy.CellVariable: Updated solution variable.
	"""
	
	#Get image resolution and center of geometry
	res=simulation.ICimg.shape[0]
	center=simulation.embryo.geometry.getCenter()
		
	#Define x/y coordinates of interpolation
	if 'quad' in simulation.embryo.analysis.process.keys():
		#Shift everything by center to fit with the mesh
		xInt = np.arange(center[0]+1, center[0]+res+1, 1)
		yInt = np.arange(center[1]+1, center[1]+res+1, 1)		
	else:
		xInt = np.arange(1, res+1, 1)
		yInt = np.arange(1, res+1, 1)
	
	#Generate interpolation function
	f=interp.RectBivariateSpline(xInt, yInt, simulation.ICimg, bbox=[None, None, None, None], kx=3, ky=3, s=0)
	
	#Getting cell centers
	x,y,z=simulation.mesh.mesh.getCellCenters()
	
	#Finding outer rim concentration
	if simulation.embryo.analysis.concRim==None:
		printWarning('concRim is not analyzed yet. Will use concentration outside of bleached region as approximation')
		
		#Grab offset and sidelength
		offset=simulation.embryo.offsetBleachedPx
		sidelength=simulation.embryo.sideLengthBleachedPx
		
		#Get indices outside of bleached square but inside masterROI
		indXSqu,indYSqu=pyfrp_idx_module.getSquareIdxImg(simulation.embryo.offsetBleachedPx,simulation.embryo.sideLengthBleachedPx,simulation.embryo.dataResPx)
		
		indX=pyfrp_misc_module.complValsSimple(simulation.embryo.masterROI.imgIdxX,indXSqu)
		indY=pyfrp_misc_module.complValsSimple(simulation.embryo.masterROI.imgIdxX,indYSqu)
		
		if 'quad' in simulation.embryo.analysis.process.keys():
			img=pyfrp_img_module.unflipQuad(np.flipud(simulation.ICimg))
		else:
			img=simulation.ICimg
		
		concRim=pyfrp_img_module.meanConc(img[indX,indY])
		
		print 'Approximate concRim = ', concRim
		
	else:	
		concRim=simulation.embryo.analysis.concRim
		
	#Set all values of solution variable to concRim
	phi.setValue(concRim)
	
	#Get Offset of image and check which nodes are inside image
	if 'quad' in simulation.embryo.analysis.process.keys():
		offset=[simulation.embryo.dataResPx/2,simulation.embryo.dataResPx/2]
		ins=pyfrp_idx_module.checkInsideImg(x,y,simulation.embryo.dataResPx/2,offset=offset)
	else:
		offset=[0,0]
		ins=pyfrp_idx_module.checkInsideImg(x,y,simulation.embryo.dataResPx,offset=offset)
		
	#Convert into indices
	ind=np.arange(len(x))
	ind=ind[np.where(ins)[0]]

	"""NOTE:  I think we need to match here indices inside image with the one of master ROI, so we don't apply
	values outside of masterROI (generally just background) to nodes that lie INSIDE image, but OUTSIDE of masterROI.
	"""
	
	if matchWithMaster:
		masterROI=simulation.embryo.getMasterROI()	
		ind=pyfrp_misc_module.matchVals(ind,masterROI.meshIdx)
		
	#Apply interpolation
	phi.value[ind]=f.ev(x[ind],y[ind])

	return phi

	

def mimic_imperfect_bleaching(slice_height,height,im_reg_ICs,conc_rim,rim,add_rim_from_radius,mesh,phi,sidelength,center,debug_opt):
	
	"""Mimic imperfect bleaching through cone approximation, return phi.
	
	 .. warning:: Not working in current version. Will be integrated in further versions again.
	
	"""
	
	#Making some variable names shorter
	sh=-slice_height
	h=height
	
	#Grabbing image resolution
	res=shape(im_reg_ICs)[0]
	
	#-----------------------------------------------
	#Finding outer rim concentration
	if conc_rim==None:
		
		conc_sum==0
		conc_num==0
		
		for i in range(res):
			for j in range(res):
				if add_rim_from_radius==1:
					if sqrt((i-center[0])**2+(j-center[1])**2)<=radius and sqrt((i-center[0])**2+(j-center[1])**2)>=radius-rim:
					
						conc_sum=conc_sum+im_reg_ICs[i,j]
						conc_num=conc_num+1
						
				elif add_rim_from_radius==0:
					if sqrt((i-center[0])**2+(j-center[1])**2)<=radius and sqrt((i-center[0])**2+(j-center[1])**2)>=res/2-rim:
						conc_sum=conc_sum+im_reg_ICs[i,j]
						conc_num=conc_num+1
		
		conc_rim=conc_sum/conc_num

	#-----------------------------------------------
	#Imperfect bleaching
	
	#Creating Grid
	x_grid=arange(0,res,1)
	y_grid=arange(0,res,1)
	
	x_grid, y_grid=meshgrid(x_grid,y_grid)
	grid_f=zeros(shape(x_grid))
	
	#Parameters of sigmoid function
	max_val=1
	max_min_val=1.25
	rate=0.01
	r_jump=sidelength/2.
	center_x=center[0]
	center_y=center[1]
	
	#Calculate linear equation describing how strong bleaching effect decreases as z increases
	m=(max_min_val-1)/((h-sh)/h)
	b=1-sh*(max_min_val-1)/(h-sh)
	
	#Grid of radius values
	r=sqrt((x_grid-center_x)**2+(y_grid-center_y)**2)
		
	#Debugging plots if needed
	if debug_opt==0:
		
		#Range over z-stacks 
		zs=arange(0,h,1)
		
		#Matrix of scaled ics
		ics_scaled=zeros((res,res,shape(zs)[0]))
			
		#Generate z-stack of slices
		j=0
		for z in zs:
			
			#Apply sigmoid function
			min_val=m*z/h+b
			grid_f=min_val+(max_val-min_val)/(1+exp(-rate*(r-r_jump)))
			ics_scaled[:,:,j]=grid_f*im_reg_ICs
				
			#Increasing index	
			j=j+1
		
		#Mesh points where the interpolation is supposed to be evaluated	
		points=zeros((shape(mesh.x)[0],3))
		points[:,0]=mesh.x
		points[:,1]=mesh.y
		points[:,2]=-mesh.z
		points=points.T
		
		#Start time of timer
		start_time_inter=time.clock()
		
		#Interpolation
		int_res=ndi.map_coordinates(ics_scaled,points,cval=conc_rim)
		
		#End time of timer
		print "interpolation took:", time.clock()-start_time_inter 
		print conc_rim
		phi.value=int_res
	
	#-----------------------------------------------
	#Some debugging outputs		
	elif debug_opt==1:
		
		#Different heights for debugging
		zs=[0,sh,h]
		
		ics_scaled=zeros((res,res,shape(zs)[0]))
		
		#Generate z-stack of slices
		j=0
		for z in zs:
			min_val=m*z/h+b
			grid_f=min_val+(max_val-min_val)/(1+exp(-rate*(r-r_jump)))
			ics_scaled[:,:,j]=grid_f*im_reg_ICs
			
				
			j=j+1
		
		#Plotting slices
		i=1
		j=0
		
		fig=plt.figure()
		fig.show()
		cflevels=linspace(0,im_reg_ICs.max(),25)

		for z in zs:
			
		
			fid=330+i
			i=i+1
			ax=fig.add_subplot(fid)
			curr_plt=ax.contourf(x_grid,y_grid,grid_f,levels=cflevels)
			plt.ylabel("z="+str(z))
			if i==2:
				plt.title("Scaling sigmoid")
			
			plt.colorbar(curr_plt)
			fid=330+i
			i=i+1
			ax=fig.add_subplot(fid)
			curr_plt=ax.contourf(ics_scaled[:,:,j],levels=cflevels)
			if i==3:
				plt.title("Scaled ICs")
			
			plt.colorbar(curr_plt)
			fid=330+i
			i=i+1
			ax=fig.add_subplot(fid)
			diff=ics_scaled[:,:,j]-im_reg_ICs
			curr_plt=ax.contourf(diff)
			if i==4:
				plt.title("Difference between scaled and img")
			plt.colorbar(curr_plt)
			
			j=j+1
		
		plt.draw()
		raw_input()
	
	return phi
	



		
	

