#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import os
import time
from optparse import OptionParser

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
#Defining parsing possibilities
#===========================================================================================================================================================================

usage = "usage: %prog [options]"
parser = OptionParser(usage)
     
parser.add_option("-l", "--file",action="store", type='string', dest="fn", help="Input file")
parser.add_option("-a", "--analyze", type='int',  dest="a", help="Analyze?")
parser.add_option("-s", "--simulate",type='int',  dest="s", help="Simulate?")
parser.add_option("-p", "--pin" ,type='int',  dest="p", help="Pin?")
parser.add_option("-f", "--fit",type='int',  dest="f", help="Fit?")
parser.add_option("-n", "--n_emb",type='int', dest="n", help="n_emb?")

parser.set_defaults(fn=None)
parser.set_defaults(a=1)
parser.set_defaults(s=1)
parser.set_defaults(p=1)
parser.set_defaults(f=1)
parser.set_defaults(n=0)

#Getting parameters
(options, args) = parser.parse_args()     

print options.fn
print options.a 

fn=str(options.fn)
analyze=int(options.a)
simulate=int(options.s)
pin=int(options.p)
fit=int(options.f)
n_emb=int(options.n)

print "================================================================================"
print "Will analyze dataset " + str(n_emb) + " in file " + fn
print "Will do the following:"
print "Analyze: ", analyze
print "Simulate: ", simulate
print "Pin: ", pin
print "Fit: ", fit
print "Embryos: ", n_emb
print "================================================================================"

#===========================================================================================================================================================================
#Some flags what to do
#===========================================================================================================================================================================


preproc=0
find_bound=0

#===========================================================================================================================================================================
#Parameters
#===========================================================================================================================================================================

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating embryo object

mol=molecule("blub")
mol=mol.load_molecule(fn)
emb=mol.embryos[n_emb]

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Analyze embryo

if analyze==1:
	
	emb.debug_analysis=0
	emb = pyfrp_img.analyze_dataset(emb)	
	print "Done Image Analysis"
	mol.save_molecule(fn)	
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Preproc for ICs 
	
if preproc==1:
	emb.rad_step_px=5

	emb = pyfrp_sim.apply_data_ics(emb)	
	emb.reg_mesh_opt=0
	print "Done Preproc"
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Simulate for scaling solution 

if simulate==1:

	#emb.conc_rim=None
	
	#Some settings
	#Plots
	emb.plot_surf=0
	emb.plot_cont=0
	emb.plot_conc=0
	emb.plot_wire=0
	emb.plot_all=0
	
	#Create Mesh on the fly
	emb.usemesh=0
	emb.geometry="Cylinder"
	emb.volSize_px=25
	
	#Debug?
	emb.debug_simulation=0
	
	#Add rim?
	emb.add_rim_from_radius=1
	
	#Apply ideal conditions
	emb.apply_data=0
	
	#Timesteps
	emb.steps_sim=1000
	
	#Converty to logarithmic scale
	emb.to_log_scale()
	
	#run PDE solver
	emb=pyfrp_sim.simulate_diff_react(emb)
	mol.save_molecule(fn)
	print "Done Simualtion"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Pin concentration profiles
	
if pin==1:	
	emb=pyfrp_fit.pin_conc(emb)	
	print "Done Pinning"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fit sim to data
	
if fit==1:
	
	#Deleting all fits
	emb.fits=[]
	emb.debug_fit=0
	
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
	emb.fits[0].x0=[40,0,0]
	emb.fits[0].equ_on=1
	emb.fits[0].fit_cut_off_t=0
	
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
	emb.fits[1].x0=[40,0,0]
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
	emb.fits[2].x0=[40,0,0]
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
	emb.fits[3].x0=[40,0,0]
	emb.fits[3].equ_on=0
	emb.fits[3].fit_cut_off_t=1
	emb.fits[3].cut_off_t=50

	#Perform fits
	emb=pyfrp_fit.parm_fitting(emb,0)
	emb=pyfrp_fit.parm_fitting(emb,1)
	emb=pyfrp_fit.parm_fitting(emb,2)
	emb=pyfrp_fit.parm_fitting(emb,3)
	
	print "Done Fitting"
	mol.save_molecule(fn)

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

