#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

from fipy import *
import math 
import os
from numpy import *
from pyfrp_sim_module import *
from pyfrp_misc_module import *
from pyfrp_img_module import *
from pyfrp_fit_module import *
from pyfrp_stats_module import *
from embryo import *
import time
import matplotlib.pyplot as plt


#Some flags what to do
load_emb=1
analyze=0
preproc=0
simulate=1
pin=1
fit=1

#Define dataset
dataset="FRAP_10kDa-A488-Dextran_05052014"
exp="exp9"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating embryo object
emb=embryo("blub","frap")
emb.name="10kDa_2014_05_05_in_water"
emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"

if load_emb==0:
	
	#Dataset
	emb.fn_datafolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/post"
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
	emb.fn_preimage="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/pre/10kDa-A488-Dex-05mgml-water_720um_100um_20um-depth_exp9_pre000.tif"
	
	emb.data_res_px=512.
	emb.conv_fact=517.17/emb.data_res_px
	emb.framerate=1
	emb.add_rim_from_radius=1
	emb.norm_by_pre=1
	
	emb.file_list=get_sorted_folder_list(emb.fn_datafolder,emb.data_ft)
	#emb.file_list=emb.file_list[:1]

	emb.nframes=shape(emb.file_list)[0]
	emb.tstart=0
	emb.tend=emb.framerate*(emb.nframes-1)
	emb.tvec_data=linspace(emb.tstart,emb.tend,emb.nframes)
	emb.steps_data=emb.nframes
	emb.dt_data=emb.framerate

	emb.side_length_bleached_mu=141.7/2
	emb.slice_depth_mu=80.
	emb.radius_embr_mu=305.
	emb.slice_width_mu=2.5
	emb.cylinder_radius_mu=emb.radius_embr_mu
	emb.cylinder_height_mu=100.
	emb.rim=50
	
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

	emb.img_in_domain=0
	emb.D=400.
	emb.usemesh=1
	emb.fn_mesh=emb.fn_resultfolder+"/meshfile_slice.msh"
	
	
else:
	fn_load=emb.fn_resultfolder+"/"+emb.name+".pk"
	emb=emb.load_embryo(fn_load)
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
	emb.fn_mesh=emb.fn_resultfolder+"/meshfile_slice3.msh"
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Analyze embryo

if analyze==1:
	emb.debug_analysis=0
	emb.norm_by_pre=1
	emb = analyze_dataset(emb)
	
	emb.save_embryo(None)
	print "Done Image Analysis"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Simulate for scaling solution 

if simulate==1:
			#emb.im_reg_ICs=gen_fake_IC(512,0,1,emb.offset_bleached_px,emb.side_length_bleached_px,emb.radius_embr_px,emb.center_embr_px,emb.rim,emb.add_rim_from_radius,0)
	squ_shs=[]
	out_shs=[]
	slice_shs=[]
	
	fig=plt.figure()
	fig.show()
	fig.subplots_adjust(right=0.8)
	ax_squ=fig.add_subplot(221)
	ax_out=fig.add_subplot(222)
	ax_slice=fig.add_subplot(223)
	
	cmap=plt.get_cmap("jet")
	
	shs=linspace(5,105,10)
	i=0
	for sd in shs: 
	
		emb.slice_depth_mu=sd
		emb.slice_depth_px=emb.slice_depth_mu/emb.conv_fact
		emb.slice_height_px=[-emb.slice_depth_px]
		
		#Some settings
		#emb=crop_conc_rim_from_pre(emb)
		
		#print emb.conc_rim
		#raw_input()
		emb.conc_rim=1.3
		emb.steps_sim=3000
		emb.tend=300
		emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)
		emb.apply_data=3
		emb.plot_surf=0
		emb.plot_cont=0
		emb.plot_conc=0
		emb.debug_simulation=0
		emb.add_rim_from_radius=1
		emb.D=80/(emb.conv_fact**2)
		
		#Load mesh
		emb.mesh=GmshImporter3D(emb.fn_mesh)
		
		#run PDE solver
		emb=diff_reg(emb)
		
		print "Done Simualtion with slice_height_px=", emb.slice_depth_mu
	
		squ_shs.append(list(emb.squ_av_d))
		out_shs.append(list(emb.out_av_d))
		slice_shs.append(list(emb.slice_av_d))
	
		color=cmap(float(i)/float(len(shs)))
	
		i=i+1
		
		ax_squ.plot(emb.tvec_sim,squ_shs[-1],color=color)
		ax_out.plot(emb.tvec_sim,out_shs[-1],color=color,label=round(sd))
		ax_slice.plot(emb.tvec_sim,slice_shs[-1],color=color)
		
	ax_squ.set_title("Squ")	
	ax_slice.set_title("Slice")
	ax_out.set_title("Out")
	
	ax_out.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	
	plt.draw()
	raw_input()
	