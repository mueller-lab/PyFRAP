from fipy import *
import math 
import os
from numpy import *
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_img_module as pyfrp_img
import pyfrp_fit_module as pyfrp_fit
import pyfrp_gmsh_module as pyfrp_gmsh
import pyfrp_plot_module as pyfrp_plot

from embryo import *
import time

#~~~~~~~~~~~~~~~~~~~~~~~~
#Creating embryo object
#~~~~~~~~~~~~~~~~~~~~~~~~

emb=embryo("blub","frap")

#Filepath
cwd=os.getcwd()
emb.fn_datafolder=cwd+"/testdata/"
emb.file_list=get_sorted_folder_list(emb.fn_datafolder,emb.data_ft)

#For perfomance purpose, take only every third image
emb.file_list=emb.file_list[::3]

#Dataset Settings
emb.nframes=shape(emb.file_list)[0]
emb.framerate=30
emb.tstart=0
emb.tend=emb.framerate*(emb.nframes-1)
emb.tvec_data=linspace(emb.tstart,emb.tend,emb.nframes)
emb.steps_data=emb.nframes
emb.dt_data=emb.framerate

#Simulation Settings
emb.steps_sim=500
emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)

#Some flags what to do
load_emb=1
analyze=0
simulate=1
pin=0
fit=0


if load_emb==1:
	emb=emb.load_embryo(cwd+"/testdata/"+"testembryo.pk")

if analyze==1:
	#Image Analysis
	emb = pyfrp_img.analyze_dataset(emb)
	emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")
	print 
	print "Done with data analysis"
	#raw_input("Press Enter")
	
	
if simulate==1:
	#Simulation
	
	emb.steps_sim=100
	emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)
	
	emb=pyfrp_sim.simulate_diff_react(emb)
	emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")
	#raw_input()
if pin==1:
	#Pinning
	emb=pyfrp_fit.pin_conc(emb)
	emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")

if fit==1:
	#Fitting
	emb.fits=[]
	
	#Fit equ, pinned
	emb.add_fit(0,"fit equ pinned","default")

	emb.fits[0].opt_meth="Nelder-Mead"
	emb.fits[0].fit_prod=1
	emb.fits[0].fit_degr=0
	emb.fits[0].fit_pinned=1
	emb.fits[0].fit_slice=1
	emb.fits[0].fit_squ=1
	emb.fits[0].fit_out=0
	emb.fits[0].UB_D=400
	emb.fits[0].x0=[8,7.2,5.37]
	emb.fits[0].equ_on=1
	
	
	emb.add_fit(1,"fit unpinned","default")			

	emb.fits[1].opt_meth="Nelder-Mead"
	emb.fits[1].fit_prod=1
	emb.fits[1].fit_degr=0
	emb.fits[1].fit_slice=1
	emb.fits[1].fit_squ=0
	emb.fits[1].fit_out=1
	emb.fits[1].fit_pinned=0
	emb.fits[1].UB_D=400
	emb.fits[1].x0=[10,3.87,5.37]
	emb.fits[1].equ_on=1
	
	emb.debug_fit=0
	
	emb=pyfrp_fit.parm_fitting(emb,0)
	emb.debug_fit=0
	
	emb=pyfrp_fit.parm_fitting(emb,1)
	
	emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")

	#Printing final results	
	emb.fits[0].print_results()
	emb.fits[0].plot_fit_pinned()

	emb.fits[1].print_results()
	emb.fits[1].plot_fit_unpinned()
	