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

	
#-------------------------------------------------------------------------------------------------------------------------------
#Settings
#-------------------------------------------------------------------------------------------------------------------------------

#Some flags what to do
sim_movie=1
fit_movie=0
ideal=1
fish=0


#Choosing dataset

emb=embryo("blub","default")

if ideal==1:
	fn_load="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Ben/1um-A488-10kDa_1-and-a-half-percent-agarose_exp1/results/10kDa_1-and-a-half-percent-agarose_exp1.pk"
	emb=emb.load_embryo(fn_load)
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Ben/1um-A488-10kDa_1-and-a-half-percent-agarose_exp1/results/"
	emb.fn_mesh=emb.fn_resultfolder+"meshfile_field.msh"
	
	print emb.slice_height_px
	
if fish==1:	
	fn_load="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/myembryo.pk"
	emb=emb.load_embryo(fn_load)
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/"
	emb.fn_mesh=emb.fn_resultfolder+"meshfile.msh"


if sim_movie==1:
	emb.mesh=GmshImporter3D(emb.fn_mesh)
	#emb.slice_width_px=20.
	emb.img_in_domain=1
	emb.add_rim_from_radius=1
	emb.conc_rim=1.0
	emb.steps_sim=1000
	emb.D=20.41
	emb.debug_simulation=0
	emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)

	emb.plot_cont=1
	emb.plot_surf=0
	emb.debug_reg=0
	emb.reg_mesh_opt=0
	emb.res_reg=100
	
	
	

	emb=diff_reg(emb)
	
if fit_movie==1:

	#Deleting all fits
	emb.fits=[]
	emb.debug_fit=1
	
	#Fit equ, pinned
	emb.add_fit(0,"default")

	emb.fits[0].opt_meth="Nelder-Mead"
	emb.fits[0].fit_prod=0
	emb.fits[0].fit_degr=0
	emb.fits[0].fit_pinned=0
	emb.fits[0].fit_slice=1
	emb.fits[0].fit_squ=1
	emb.fits[0].fit_out=0
	emb.fits[0].UB_D=400
	emb.fits[0].x0=[50,0,0]
	emb.fits[0].equ_on=0	
	
	#emb.fits[0].x0[2]=5.37
	
	
	emb=dim_fitting(emb,0)
	
	
	
	


