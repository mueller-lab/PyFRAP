#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import os
import time
from optparse import OptionParser
import copy as cpy

#Numpy/Scipy
from numpy import *

#PyFRAP Modules
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_img_module as pyfrp_img
import pyfrp_fit_module as pyfrp_fit
import pyfrp_gmsh_module as pyfrp_gmsh
import pyfrp_plot_module as pyfrp_plot
from embryo import *
from molecule import *

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

fn_mol="../Ideal_FRAP/Gary/testmol.pk"

sim=False

if sim:

	#Load some molecule object with some data in there
	mol=molecule("blub")
	mol=mol.load_molecule(fn_mol)

	#Grab some embryo
	emb=mol.embryos[2]

	#Create new molecule
	mol2=molecule("testmol2")
 
	#Add embryo twice to newly created molecule
	mol2.embryos.append(emb)
	embcpy=cpy.deepcopy(emb)
	mol2.embryos.append(embcpy)

	emb0=mol2.embryos[0]
	emb1=mol2.embryos[1]

	#Simulate once with D=10 and 300 steps
	emb0.D=0.5
	emb0.steps_sim=300
	emb0.to_lin_scale()
	emb0=pyfrp_sim.simulate_diff_react(emb0)

	#Simulate once with D=10 and 3000 steps
	emb1.D=10.
	emb1.steps_sim=3000
	
	emb1.tvec_sim=linspace(emb1.tstart,3000.,emb1.steps_sim)
	print min(emb1.tvec_sim), max(emb1.tvec_sim)
	emb1=pyfrp_fit.get_opt_tvec_sim_px(emb1,100)
	print min(emb1.tvec_sim), max(emb1.tvec_sim)
	
	#raw_input()
	
	#emb1.to_lin_scale()
	emb1=pyfrp_sim.simulate_diff_react(emb1)

	#Now assume that data of emb0 is the actually simulated data
	emb0.squ_av_data_d=list(emb0.squ_av_d)
	emb0.out_av_data_d=list(emb0.out_av_d)
	emb0.slice_av_data_d=list(emb0.slice_av_d)

	#And put the simulated data of emb1 as the simulated data of emb0
	emb0.squ_av_d=list(emb1.squ_av_d)
	emb0.out_av_d=list(emb1.out_av_d)
	emb0.slice_av_d=list(emb1.slice_av_d)
	emb0.tvec_sim=emb1.tvec_sim.copy()
	
	emb0.D=10.
	
	#Save molecule
	mol2.save_molecule("../Ideal_FRAP/Gary/testmol2.pk")

	#Plot and see if everything worked
	#emb0.plot_sim_data()
	
else:
	#Load new molecule
	mol2=molecule("testmol2")
	mol2=mol2.load_molecule("../Ideal_FRAP/Gary/testmol2.pk")
	
	emb0=mol2.embryos[0]
	emb0.D=10.
	#Plot and see if everything worked
	emb0.plot_sim_data()
	
#Pin emb0	
emb0.debug_pinning=1
emb0=pyfrp_fit.pin_conc(emb0)	

emb0.plot_pinned()

#Deleting all fits
emb0.fits=[]
emb0.debug_fit=1

#Fit equ, pinned
emb0.add_fit(0,"")

emb0.fits[0].opt_meth="Constrained Nelder-Mead"
emb0.fits[0].fit_prod=0
emb0.fits[0].fit_degr=0
emb0.fits[0].fit_pinned=1
emb0.fits[0].fit_slice=1
emb0.fits[0].fit_squ=1
emb0.fits[0].fit_out=0
emb0.fits[0].UB_D=400.
emb0.fits[0].LB_D=0.01
emb0.fits[0].x0=[10,0,0]
emb0.fits[0].equ_on=1
emb0.fits[0].fit_cut_off_t=0
emb0.fits[0].opt_tol=0.0000000001

#Fit equ, pinned TNC
emb0.add_fit(1,"")

emb0.fits[1].opt_meth="TNC"
emb0.fits[1].fit_prod=0
emb0.fits[1].fit_degr=0
emb0.fits[1].fit_pinned=1
emb0.fits[1].fit_slice=1
emb0.fits[1].fit_squ=1
emb0.fits[1].fit_out=0
emb0.fits[1].UB_D=400.
emb0.fits[1].LB_D=0.01
emb0.fits[1].x0=[20,0,0]
emb0.fits[1].equ_on=1
emb0.fits[1].fit_cut_off_t=0

emb0=pyfrp_fit.parm_fitting(emb0,0)
emb0=pyfrp_fit.parm_fitting(emb0,1)

print "-----------------------"
print "pinned eqon CNM:"
emb0.fits[0].print_results()
print "D_opt_px=", emb0.fits[0].D_opt_px
emb0.fits[0].plot_fit_pinned()
print "-----------------------"
print "pinned eqon TNC:"
emb0.fits[1].print_results()
print "D_opt_px=", emb0.fits[1].D_opt_px
emb0.fits[1].plot_fit_pinned()

