#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import os
import time

#Numpy/Scipy
from numpy import *

#PyFRAP Modules
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_img_module as pyfrp_img
import pyfrp_fit_module as pyfrp_fit
import pyfrp_gmsh_module as pyfrp_gmsh
from embryo import *

#===========================================================================================================================================================================
#Some flags what to do
#===========================================================================================================================================================================

load_emb=0
analyze=1
preproc=0

find_bound=0

simulate=0
pin=0
fit=0

#===========================================================================================================================================================================
#Parameters
#===========================================================================================================================================================================

#Define dataset
#dataset=""
#exp="exp10"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating embryo object
emb=embryo("blub","frap")
emb.name=dataset+"_"+exp
emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"

if load_emb==0:
	
	#Dataset
	emb.fn_datafolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/post"
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
	
	
	if find_bound==1:
		emb.norm_by_pre=0
	else:
		emb.norm_by_pre=1
	
	if emb.norm_by_pre==1:
		prefolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/pre/"
		prelist=pyfrp_misc.get_sorted_folder_list(prefolder,emb.data_ft)
		emb.fn_preimage=prefolder+prelist[0]
		
	emb.data_res_px=512.

	emb.conv_fact=485.82/emb.data_res_px
	
	
	emb.framerate=1

	
	emb.file_list=pyfrp_misc.get_sorted_folder_list(emb.fn_datafolder,emb.data_ft)
	
	if find_bound==1:
		emb.file_list=emb.file_list[:1]

	emb.nframes=shape(emb.file_list)[0]
	emb.tstart=0
	emb.tend=emb.framerate*(emb.nframes-1)
	emb.tvec_data=linspace(emb.tstart,emb.tend,emb.nframes)
	emb.steps_data=emb.nframes
	emb.dt_data=emb.framerate

	emb.side_length_bleached_mu=141.7/2
	emb.slice_depth_mu=30.
	emb.radius_embr_mu=290.
	emb.slice_width_mu=2.5
	emb.cylinder_radius_mu=emb.radius_embr_mu
	emb.cylinder_height_mu=100.
	emb.rim=50
	
	emb.gaussian=1
	emb.gaussian_sigma=3.
	
	emb.debug_analysis=1
	
	#Geometry
	emb.geometry="Cylinder"
	emb.slice_width_px=emb.slice_width_mu/emb.conv_fact
	emb.slice_depth_px=emb.slice_depth_mu/emb.conv_fact
	emb.slice_height_px=[-emb.slice_depth_px]
	emb.slice_bottom=0
	emb.center_embr_px=[emb.data_res_px/2,emb.data_res_px/2]
	

	
	emb.side_length_bleached_px=2*emb.side_length_bleached_mu/emb.conv_fact
	emb.offset_bleached_px=[emb.center_embr_px[0]-emb.side_length_bleached_px/2, emb.center_embr_px[1]-emb.side_length_bleached_px/2]
	emb.radius_embr_px=emb.radius_embr_mu/emb.conv_fact
	emb.cylinder_height_px=emb.cylinder_height_mu/emb.conv_fact
	emb.cylinder_radius_px=emb.cylinder_radius_mu/emb.conv_fact

	#Preproc
	emb.rad_step_px=100
	
	#Simulation 
	emb.steps_sim=3000
	emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)
	emb.debug_simulation=0
	emb.avg_mode=0
	emb.add_rim_sim=0
	emb.avg_outer=0
	emb.avg_inner=0
	emb.avg_all=0
	emb.avg_pocket=0
	emb.avg_small=0
	emb.int_steps=400.
	emb.conc_rim=1.0
	emb.apply_data=3
	emb.volSize_px=11
	emb.D=400/(emb.conv_fact**2)
	
	
	emb.img_in_domain=0
	
	emb.add_rim_from_radius=1
	emb.usemesh=1
	emb.fn_mesh=emb.fn_resultfolder+"/meshfile.msh"
	
	
	emb.fn_geo=emb.fn_resultfolder+"/meshfile.geo"
	
	#if find_bound==0:
	
		#pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"radius",emb.radius_embr_px)
		#pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"height",emb.cylinder_height_px)
		#pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"center_x",emb.center_embr_px[0])
		#pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"center_y",emb.center_embr_px[1])
		#pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"volSize_px",emb.volSize_px)
		#os.system("gmsh -3 " + emb.fn_geo)
		
		#Saving embryo
		#emb.save_embryo(None)
		
		
		
		
		#print "Generated Embryo Object."
		
		#raw_input()

else:

	#fn_load=emb.fn_resultfolder+"/"+emb.name+".pk"
	fn_load="Y:\People\Gary\03 - FRAP + FCS\20150401\beads_embryo.pk"
	emb=emb.load_embryo(fn_load)
	#emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
	
	#emb.fn_mesh=emb.fn_resultfolder+"/meshfile.msh"
	

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Analyze embryo

if analyze==1:
	if find_bound==1:
		emb.debug_analysis=1
	else:
		emb.debug_analysis=0
		
	emb.debug_analysis=1
	emb = pyfrp_img.analyze_dataset(emb)
	
	#emb.save_embryo(None)
	print "Done Image Analysis"
	raw_input()
	#raw_input()
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Preproc for ICs 
	
if preproc==1:
	emb.rad_step_px=5
	

	emb = pyfrp_sim.apply_data_ics(emb)
	
	emb.reg_mesh_opt=0
	
	#emb.save_embryo(None)
	
	print "Done Preproc"
	raw_input()
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Simulate for scaling solution 

if simulate==1:
	
	#Some settings
	emb.conc_rim=1.0
	emb.apply_data=3
	emb.plot_surf=0
	emb.plot_cont=0
	emb.plot_conc=0
	emb.plot_wire=0
	emb.plot_all=0
	
	emb.usemesh=0
	emb.geometry="Fish"
	
	emb.debug_simulation=0
	emb.add_rim_from_radius=1
	
	emb.steps_sim=400
	emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)
	
	#run PDE solver
	emb=pyfrp_sim.simulate_diff_react(emb)
	
	#emb.save_embryo(None)
	print "Done Simualtion"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pin concentration profiles
	
if pin==1:
	
	
	emb=pyfrp_fit.pin_conc(emb)	
	#emb.save_embryo(None)
	print "Done Pinning"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fit sim to data
	
if fit==1:
	
	#Deleting all fits
	emb.fits=[]
	emb.debug_fit=1
	
	#Fit equ, pinned
	emb.add_fit(0,"","default")

	emb.fits[0].opt_meth="Nelder-Mead"
	emb.fits[0].fit_prod=0
	emb.fits[0].fit_degr=0
	emb.fits[0].fit_pinned=1
	emb.fits[0].fit_slice=0
	emb.fits[0].fit_squ=1
	emb.fits[0].fit_out=0
	emb.fits[0].UB_D=400
	emb.fits[0].x0=[50,0,0]
	emb.fits[0].equ_on=1
	
	#Fit unpinned
	emb.add_fit(1,"","default")			

	emb.fits[1].opt_meth="Nelder-Mead"
	emb.fits[1].fit_prod=0
	emb.fits[1].fit_degr=0
	emb.fits[1].fit_slice=0
	emb.fits[1].fit_squ=1
	emb.fits[1].fit_out=0
	emb.fits[1].fit_pinned=0
	emb.fits[1].UB_D=400
	emb.fits[1].x0=[82,0,0]
	emb.fits[1].equ_on=0
	
	#Fit equ, pinned, cut_off
	emb.add_fit(2,"","default")			

	emb.fits[2].opt_meth="Nelder-Mead"
	emb.fits[2].fit_prod=0
	emb.fits[2].fit_degr=0
	emb.fits[2].fit_pinned=1
	emb.fits[2].fit_slice=0
	emb.fits[2].fit_squ=1
	emb.fits[2].fit_out=0
	emb.fits[2].UB_D=400
	emb.fits[2].x0=[50,0,0]
	emb.fits[2].equ_on=1
	emb.fits[2].fit_cut_off_t=1
	emb.fits[2].cut_off_t=50
	
	#Fit cut_off
	emb.add_fit(3,"","default")			

	emb.fits[3].opt_meth="Nelder-Mead"
	emb.fits[3].fit_prod=0
	emb.fits[3].fit_degr=0
	emb.fits[3].fit_slice=0
	emb.fits[3].fit_squ=1
	emb.fits[3].fit_out=0
	emb.fits[3].fit_pinned=0
	emb.fits[3].UB_D=400
	emb.fits[3].x0=[50,0,0]
	emb.fits[3].equ_on=0
	emb.fits[3].fit_cut_off_t=1
	emb.fits[3].cut_off_t=50
	
	print emb.steps_data
	
	#Perform fits
	emb=pyfrp_fit.parm_fitting(emb,0)
	emb=pyfrp_fit.parm_fitting(emb,1)
	emb=pyfrp_fit.parm_fitting(emb,2)
	emb=pyfrp_fit.parm_fitting(emb,3)
	
	#emb.save_embryo(None)

	print "Done Fitting"
	raw_input()

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

