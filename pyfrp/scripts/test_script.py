#! /usr/bin/python

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

test="sim2"
print "test=",test


if test=="sim2":
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#Settings
	#-------------------------------------------------------------------------------------------------------------------------------
	
	#Some flags what to do
	
	save_scal=0
	save_real=0
	fit_sol=0
	pin_sol=0
	
	load_pm=0
	
	
	#Choosing dataset
	#emb2=embryo("otherembryo","default")
	emb=embryo("blub","frap")
	
	
	fn_load="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/myembryo.pk"
	emb=emb.load_embryo(fn_load)
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/"
	
	#emb.plot_sim_data()
	#emb.fits[0].plot_fit_pinned()
	#raw_input()
	
	
	
	#emb2.copy_embryo(emb)
	
	
	#for i in range(emb.fit_number-1):
		#print emb.fits[i].opt_meth
		#print emb.fits[i].success
		#print emb.fits[i].D_opt_mu, emb.fits[i].prod_opt_scaled, emb.fits[i].degr_opt_scaled 
		
		#emb.fits[i].plot_fit()
	
	#emb.print_embryo()
	#emb.fits[0].print_fit()
	#emb.plot_sim_data()
	#emb.plot_pinned()
	
	emb.mode="frap"
	
	
	#Load mesh
	#emb.mesh=GmshImporter3D(emb.fn_mesh)
		
	#-------------------------------------------------------------------------------------------------------------------------------
	#Analyze dataset
	#-------------------------------------------------------------------------------------------------------------------------------
	
	if save_real==1:
		#Save data
		emb=analyze_dataset(emb)
		
		emb.save_embryo()
		print "done analysis, press enter"
		raw_input()
	
		
		
	if load_pm==1:
		#Load data and save to embryo file
		
		fn_save=fn_results+"/squ_av_pm_d.txt"
		emb.squ_av_data_pm_d=loadtxt(fn_save)	
		
		fn_save=fn_results+"/out_av_pm_d.txt"
		emb.out_av_data_pm_d=loadtxt(fn_save)
		
		fn_save=fn_results+"/slice_av_pm_d.txt"
		emb.slice_av_data_pm_d=loadtxt(fn_save)
		
		fn_save=fn_results+"/all_av_pm_d.txt"
		emb.all_av_data_pm_d=loadtxt(fn_save)
	
		fn_save=fn_results+"/inner_av_pm_d.txt"
		emb.inner_av_data_pm_d=loadtxt(fn_save)
	
		fn_save=fn_results+"/outer_av_pm_d.txt"
		emb.outer_av_data_pm_d=loadtxt(fn_save)
	
		fn_save=fn_results+"/out_pocket_av_pm_d.txt"
		emb.out_pocket_av_data_pm_d=loadtxt(fn_save)
		
		emb.save_embryo()
		raw_input()
		


	#-------------------------------------------------------------------------------------------------------------------------------
	#Preporcessing for ICs
	#-------------------------------------------------------------------------------------------------------------------------------
	
	#if save_scal==1:
		#emb = apply_data_ics(emb)
		#emb.save_embryo()
		#print "done preprocessing, press enter"
		#raw_input()
		#raw_input()
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#Creating scaling solution or fake dataset
	#-------------------------------------------------------------------------------------------------------------------------------

	#if save_scal==1:
		#emb.plot_cont=1
		#emb.plot_surf=1
		#emb.debug_reg=1
		#emb.reg_mesh_opt=0
		#emb.res_reg=100
		
		#emb.debug_simulation=1
		#emb=diff_reg(emb)
		#print "done simulation, press enter"
		#raw_input()
		##emb.save_embryo()
		
	#-------------------------------------------------------------------------------------------------------------------------------
	#Pin solution between 0 and 1
	#-------------------------------------------------------------------------------------------------------------------------------
	
	if pin_sol==1:
	
		emb.debug_pinning=1
		emb=pin_conc(emb)
	
		emb.save_embryo()
	
	#-------------------------------------------------------------------------------------------------------------------------------
	#Fit simulation to data
	#-------------------------------------------------------------------------------------------------------------------------------
	 
	#Fitting
	if fit_sol==1:
		emb.fits=[]
		emb.fit_number=0
		methds=["Nelder-Mead"]
		for meth in methds:
			emb.add_fit(emb.fit_number,"pinned e","default")
			
			emb.fits[emb.fit_number-1].opt_meth=meth
		
			emb.debug_fit=0
			emb.fits[emb.fit_number-1].x0[0]=50
			emb.fits[emb.fit_number-1].x0[2]=5.37

			emb=dim_fitting(emb,emb.fit_number-1)
	
			print "Done with testing", meth
	
	
	#emb.fits[0].plot_fit_pinned()
	#emb.fits[0].print_results()
	##emb.save_embryo(None)
	
	emb.name="emb4"
	emb.fits=[]
	emb.fit_number=0
	emb.debug_pinning=0
	fn_save="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/goodfit.pk"
	
	emb.save_embryo(fn_save)
	
	
	
	
	
	
