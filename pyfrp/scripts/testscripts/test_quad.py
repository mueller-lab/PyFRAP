#===========================================================================================================================================================================
#Script to test quadrant simulation
#===========================================================================================================================================================================

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import os
import time
from optparse import OptionParser

#Numpy/Scipy
from numpy import *
from scipy.interpolate import griddata

#PyFRAP Modules
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_img_module as pyfrp_img
import pyfrp_fit_module as pyfrp_fit
import pyfrp_gmsh_module as pyfrp_gmsh
import pyfrp_plot_module as pyfrp_plot
import pyfrp_integration_module as pyfrp_int
from embryo import *
from molecule import *

#Matplotlib
import matplotlib.pyplot as plt

#Image Processing
import skimage.filters.rank as skirank
import skimage.morphology as skimorph


#===========================================================================================================================================================================
#Script
#===========================================================================================================================================================================

analyze=True
simulate=True
fit=True
pin=True

fn_mol="../Ideal_FRAP/Gary/20150812/Analysis_beads.pk"

fn_out="../Ideal_FRAP/Gary/20150812/Analysis_beads_norm.pk"

#Start timer
start_time=time.clock()

#----------------------------
#Load molecule
mol=molecule("bla")
mol=mol.load_molecule(fn_mol)

#----------------------------
#Grab first embryo
emb=mol.embryos[1]

emb.debug_analysis=False
emb.quad_red=True
emb.flip_before_process=True
emb.norm_by_pre=True
emb.gaussian=False

#emb.conv_fact=322./emb.data_res_px

emb.conv_fact=566.79/emb.data_res_px

emb.fn_preimage="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Gary/20150812/02/pre/Beads_Fscin70kDa_500nM_02_pre_t000.tif"

#----------------------------
#Analyze
if analyze:

	emb=pyfrp_img.analyze_dataset(emb)
	mol.save_molecule(fn_out)

#----------------------------
#Simulate
if simulate:

	#Plot?
	emb.plot_surf=0
	emb.plot_cont=0
	emb.plot_conc=0
	emb.plot_wire=0
	emb.plot_all=0

	#Set embryo height
	emb.cylinder_height_px=100/emb.conv_fact
	emb.slice_depth_mu=40.
	emb.slice_depth_px=emb.slice_depth_mu/emb.conv_fact
	emb.slice_height_px=[-emb.slice_depth_px]
	
	#Create Mesh on the fly
	emb.usemesh=1
	emb.geometry="Cylinder"
	emb.volSize_px=23

	#Make sure to average everything for control
	emb.avg_all=1

	pyfrp_gmsh.update_cylinder_msh("quad_cyl.geo",emb.radius_embr_px,emb.cylinder_height_px,emb.center_embr_px,emb.volSize_px)
	pyfrp_gmsh.refine_msh("quad_cyl.msh")
	#pyfrp_gmsh.refine_msh("quad_cyl.msh")

	emb.fn_mesh="quad_cyl.msh"

	#Debug?
	emb.debug_simulation=1

	#Add rim?
	emb.add_rim_from_radius=1

	#Apply ICs via interpolation
	emb.apply_data=3

	emb.D=100

	#Timesteps
	emb.steps_sim=1000

	#Reduction to quadrant
	embryo.quad_red=True


	#Converty to logarithmic scale
	emb.to_log_scale()
	#emb.to_lin_scale()

	#run PDE solver
	emb=pyfrp_sim.simulate_diff_react(emb)
	#emb.plot_sim_data()

	mol.save_molecule(fn_out)

if pin:	
	emb=pyfrp_fit.pin_conc(emb)	
	print "Done Pinning"
	
if fit:
				
	#Deleting all fits
	emb.fits=[]
	emb.debug_fit=0
	
	#Fit equ, pinned
	emb.add_fit(0,"")
	
	emb.fits[0].opt_meth="TNC"
	emb.fits[0].fit_prod=0
	emb.fits[0].fit_degr=0
	emb.fits[0].fit_pinned=1
	emb.fits[0].fit_slice=0
	emb.fits[0].fit_squ=1
	emb.fits[0].fit_out=0
	emb.fits[0].UB_D=400
	emb.fits[0].LB_D=0
	emb.fits[0].x0=[20,0,0]
	emb.fits[0].equ_on=1
	emb.fits[0].fit_cut_off_t=0
	
	#Fit unpinned
	emb.add_fit(1,"")			

	emb.fits[1].opt_meth="TNC"
	emb.fits[1].fit_prod=0
	emb.fits[1].fit_degr=0
	emb.fits[1].fit_slice=0
	emb.fits[1].fit_squ=1
	emb.fits[1].fit_out=0
	emb.fits[1].fit_pinned=0
	emb.fits[1].UB_D=400
	emb.fits[0].LB_D=0
	emb.fits[1].x0=[20,0,0]
	emb.fits[1].equ_on=0
	
	#Fit equ, pinned, cut_off
	emb.add_fit(2,"")			

	emb.fits[2].opt_meth="TNC"
	emb.fits[2].fit_prod=0
	emb.fits[2].fit_degr=0
	emb.fits[2].fit_pinned=1
	emb.fits[2].fit_slice=0
	emb.fits[2].fit_squ=1
	emb.fits[2].fit_out=0
	emb.fits[2].UB_D=400
	emb.fits[0].LB_D=0
	emb.fits[2].x0=[20,0,0]
	emb.fits[2].equ_on=1
	emb.fits[2].fit_cut_off_t=1
	emb.fits[2].cut_off_t=50
	
	#Fit cut_off
	emb.add_fit(3,"")			

	emb.fits[3].opt_meth="TNC"
	emb.fits[3].fit_prod=0
	emb.fits[3].fit_degr=0
	emb.fits[3].fit_slice=0
	emb.fits[3].fit_squ=1
	emb.fits[3].fit_out=0
	emb.fits[3].fit_pinned=0
	emb.fits[3].UB_D=400
	emb.fits[0].LB_D=0
	emb.fits[3].x0=[20,0,0]
	emb.fits[3].equ_on=0
	emb.fits[3].fit_cut_off_t=1
	emb.fits[3].cut_off_t=50

	#Perform fits
	emb=pyfrp_fit.parm_fitting(emb,0)
	emb=pyfrp_fit.parm_fitting(emb,1)
	emb=pyfrp_fit.parm_fitting(emb,2)
	emb=pyfrp_fit.parm_fitting(emb,3)
	
	print "Done Fitting"
	mol.save_molecule(fn_out)

print "Done after:", time.clock()-start_time


print "-----------------------"
print "pinned eqon:"
emb.fits[0].print_results()
emb.fits[0].plot_fit_pinned()
print "-----------------------"
print "unpinned eqoff:"
emb.fits[1].print_results()
emb.fits[1].plot_fit_unpinned()
print "-----------------------"
print "pinned eqon cut_off:"
emb.fits[2].print_results()
emb.fits[2].plot_fit_pinned()
print "unpinned eqoff cut_off:"
emb.fits[3].print_results()
emb.fits[3].plot_fit_unpinned()	
