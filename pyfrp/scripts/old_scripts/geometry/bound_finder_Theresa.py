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
import pyfrp_plot_module as pyfrp_plot
import pyfrp_subwin

from embryo import *

#===========================================================================================================================================================================
#Some flags what to do
#===========================================================================================================================================================================

load_emb=0
gmsh=False

#===========================================================================================================================================================================
#Parameters
#===========================================================================================================================================================================

#Define dataset
dataset="FITC-dextran-40kDa_invitrogen_1uM_H2O_20140801"
exp="exp10"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating embryo object
emb=embryo("blub","frap")
emb.name=dataset+"_"+exp
emb.fn_resultfolder="/Volumes/ag-mueller/Theresa/"+dataset+"/"+exp+"/results"


if load_emb==0:
	
	#Dataset
	emb.fn_datafolder="/Volumes/ag-mueller/Theresa/"+dataset+"/"+exp+"/post"
	emb.fn_resultfolder="/Volumes/ag-mueller/Theresa/"+dataset+"/"+exp+"/results"
	
	#Turn off norm by pre so it is easier to see boundaries
	emb.norm_by_pre=0
	
	#Analysis
	emb.gaussian=0
	emb.gaussian_sigma=2
		
	emb.data_res_px=512.
	emb.conv_fact=485.82/emb.data_res_px
	
	emb.framerate=1
	
	emb.file_list=pyfrp_misc.get_sorted_folder_list(emb.fn_datafolder,emb.data_ft)
	
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
	
	emb.debug_analysis=0
	
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

	
else:

	fn_load=emb.fn_resultfolder+"/"+emb.name+".pk"
	emb=emb.load_embryo(fn_load)
	emb.fn_resultfolder="/Volumes/ag-mueller/Theresa/"+dataset+"/"+exp+"/results"
	
	emb.fn_mesh=emb.fn_resultfolder+"/meshfile.msh"
	emb.fn_geo=emb.fn_resultfolder+"/meshfile.geo"


#Call Domain selector	
selector=pyfrp_plot.boundary_selector(emb)
emb=selector.get_embryo()

#Update in Gmsh file
if gmsh:
	pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"radius",emb.radius_embr_px)
	pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"height",emb.cylinder_height_px)
	pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"center_x",emb.center_embr_px[0])
	pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"center_y",emb.center_embr_px[1])
	pyfrp_gmsh.update_parm_in_file(emb.fn_geo,"volSize_px",emb.volSize_px)
	os.system("gmsh -3 " + emb.fn_geo)

#Turn on norm_by_pre after boundary is found
emb.norm_by_pre=1
if emb.norm_by_pre==1:
	prefolder="/Volumes/ag-mueller/Theresa/"+dataset+"/"+exp+"/pre/"
	prelist=pyfrp_misc.get_sorted_folder_list(prefolder,emb.data_ft)
	emb.fn_preimage=prefolder+prelist[0]

print "Generated Embryo Object."
raw_input()

#Saving embryo
emb.save_embryo(None)	