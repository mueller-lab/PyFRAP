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
parser.add_option("-a", "--analyze", type='int',  dest="a", help="Analyze? 1=Yes, 0=No")
parser.add_option("-s", "--simulate",type='int',  dest="s", help="Simulate? 1=Yes, 0=No")
parser.add_option("-p", "--pin" ,type='int',  dest="p", help="Pin? 1=Yes, 0=No")
parser.add_option("-f", "--fit",type='int',  dest="f", help="Fit? 1=Yes, 0=No")
parser.add_option("-o", "--out",type='int',  dest="o", help="Out?")
parser.add_option("-e", "--embs",type='str',  dest="e", help="List of embryos?")
parser.add_option("-S", "--stats",type='int',  dest="S", help="Print final statistics? 1=Yes, 0=No")
parser.add_option("-O", "--fout",type='str',  dest="O", help="List of outfiles")
parser.add_option("-n", "--norm",type='int',  dest="n", help="Norm by preimage? 1=Yes, 0=No (default), 2=Keep previous setting")
parser.add_option("-g", "--geo",type='int',  dest="g", help="Geometry? 1=zebrafish, 0=cylinder (default)")
parser.add_option("-b", "--bound",type='int',  dest="b", help="Select Boundary?  1=Yes, 0=No (default)")
parser.add_option("-F", "--fix",type='int',  dest="F", help="Fix Filepaths?  1=Yes, 0=No (default)")
parser.add_option("-q", "--quad",type='int',  dest="q", help="Quadrant reduction?  1=Yes, 0=No (default)")
parser.add_option("-v", "--vol",type='int',  dest="v", help="Mesh cell volume?  volSize_px=23 (default)")
parser.add_option("-t", "--steps",type='int',  dest="t", help="Simulation steps?  steps_sim=3000 (default)")
parser.add_option("-r", "--regions",type='str',  dest="r", help="Regions to fit?  [1,0,1] (squ,out,slice) (default)")
parser.add_option("-G", "--refined",type='str',  dest="G", help="Refine slice? 1=Yes, 0=No (default)")


parser.set_defaults(fn=None)
parser.set_defaults(a=1)
parser.set_defaults(s=1)
parser.set_defaults(p=1)
parser.set_defaults(f=1)
parser.set_defaults(o=0)
parser.set_defaults(S=0)
parser.set_defaults(n=0)
parser.set_defaults(g=0)
parser.set_defaults(F=0)
parser.set_defaults(b=0)
parser.set_defaults(q=0)
parser.set_defaults(v=23)
parser.set_defaults(t=3000)
parser.set_defaults(O=None)
parser.set_defaults(r="[1,0,1]")
parser.set_defaults(G=0)

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
analyze=int(options.a)
simulate=int(options.s)
pin=int(options.p)
fit=int(options.f)
out=int(options.o)
stats=int(options.S)
embs=str(options.e)
fn_out=str(options.O)
norm_by_pre=int(options.n)
geometry=int(options.g)
bound=int(options.b)
fix=int(options.F)
quad_red=int(options.q)
steps_sim=int(options.t)
volSize_px=int(options.v)
fit_regions=str(options.r)
refine_slice=int(options.G)

print "================================================================================"
print "Will analyze all datasets specified in file " + fn
print "Will do the following:"
print "Analyze: ", analyze
print "Simulate: ", simulate
print "Pin: ", pin
print "Fit: ", fit
print "Out: ", out
print "Stats: ", stats
print "Embryos: ", embs
print "fn_out: ", fn_out
print "norm_by_pre: ", norm_by_pre
print "geometry: ", geometry
print "Fix filepaths: ", fix
print "Select boundary: ", bound
print "Quad_red: ", quad_red
print "volSize_px: ", volSize_px
print "steps_sim: ", steps_sim
print "(fit_squ,fit_out,fit_slice): ", fit_regions
print "refine_slice: ", refine_slice
print "================================================================================"

#===========================================================================================================================================================================
#Some flags what to do
#===========================================================================================================================================================================

preproc=0
find_bound=0

#===========================================================================================================================================================================
#Read in textfile with input filenames
#===========================================================================================================================================================================
fn_test="test_backup.txt"
os.system("cp " + fn + " " + fn_test)

f = open(fn, "r")
mols=[]
for line in f:
	mols.append(line.strip('\n'))
f.close()

#===========================================================================================================================================================================
#Parameters
#===========================================================================================================================================================================

#Read in embryos to analyze
if embs=='None':
	embs=len(mols)*"[],"
	embs="["+embs[:-1]+"]"
	
embs,dump=pyfrp_misc.str2list(embs)

#Read in out filenames
if fn_out=='None':
	fn_out=len(mols)*"'',"
	fn_out="["+fn_out[:-1]+"]"
	
fn_out,dump=pyfrp_misc.str2list(fn_out,dtype='str')

if len(mols)>len(fn_out):
	for i in range(len(mols)-len(fn_out)):
		fn_out.append("")
		
#Read in regions to fit
fit_regions,dump=pyfrp_misc.str2list(fit_regions,dtype='int')
fit_squ=fit_regions[0]
fit_out=fit_regions[1]
fit_slice=fit_regions[2]

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

#Loop through all molecules specified
for k,fn_mol in enumerate(mols):
	
	print "======================================================================================================================================="
	print "Processing molecule " + fn_mol
	print "======================================================================================================================================="
	
	#Creating molecule object and load molecule file
	mol=molecule("blub")
	mol=mol.load_molecule(fn_mol)
	
	#Select fn_out
	if len(fn_out[k])>3:
		fn_out_mol=fn_out[k]
	else:
		fn_out_mol=fn_mol
		
	#Empty sel_fits list
	mol.sel_fits=[]
	
	#Take all embryos if embs is empty
	if len(embs[k])==0:
		embs[k]=range(len(mol.embryos))
		
	#Adjust filepaths if necessary
	#datestr=pyfrp_misc.find_date_str(fn_mol)
	#mol=pyfrp_misc.adjust_filepaths(mol,"/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Gary/"+datestr+"/",osold="win",osnew="lin",embs=embs[k],mode='auto',debug=True)
	
	#Loop through all embryos in molecule
	for n,emb in enumerate(mol.embryos):
		
		
		#Check if we are supposed to analyze this embryo
		if n in embs[k]:
			
			print "-------------------------------------------------------------------------------------------------------------------------"
			print "Processing embryo " + emb.name + " of molecule " + fn_mol
			print "-------------------------------------------------------------------------------------------------------------------------"
			
			#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#Some settings
			
			#Analysis Settings
			emb.debug_analysis=False
			emb.quad_red=bool(quad_red)
			emb.flip_before_process=True
			emb.gaussian=False
			emb.gaussian_sigma=2.
			emb.data_offset=1
			emb.add_rim_from_radius=1
				
			#Scale rim to be 20% of embryo radius
			emb.rim=0.3*emb.radius_embr_px
			
			if norm_by_pre!=2:
				emb.norm_by_pre=bool(norm_by_pre)
			
			
			#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#Fixing filepaths
			
			#Fix preimage filepath
			if fix:
				#if not os.path.isfile(emb.fn_preimage):
					#emb.fn_preimage=pyfrp_misc.find_preimage('pre',emb.fn_datafolder,debug=False)
				#else:
					#emb.fn_preimage=pyfrp_misc.fix_fn_digits(emb.fn_preimage,3,debug=False)
				
				#Check if datafolder is there and if not try to fix it
				if os.path.isdir(emb.fn_datafolder):
					pass
				else:
					
					#Grab date
					datestr=pyfrp_misc.find_date_str(mol.name)
					
					#Find embryo number
					idxs=find_int_str(emb.name)
					nmbr=emb.name[min(idxs[-1]):]
					if len(nmbr)<2:
						nmbr="0"+nmbr
					
					#Hard adjust filepath according to Alex's file sorting
					if "nobeads" in mol.name:
						emb.fn_datafolder=os.path.split(fn_mol)[0]+"/nobeads/"+datestr+"_"+nmbr+"/recover/"
					elif "embryo" in mol.name:
						emb.fn_datafolder=os.path.split(fn_mol)[0]+"/embryo/"+datestr+"_"+nmbr+"/recover/"
					else:
						emb.fn_datafolder=os.path.split(fn_mol)[0]+"/beads/"+datestr+"_"+nmbr+"/recover/"
					
					
					emb.fn_preimage=pyfrp_misc.find_preimage('pre',emb.fn_datafolder,debug=True)	
					
					#Last checks, skip embryo if nothing worked
					if not os.path.isdir(emb.fn_datafolder):
						print "Warning, really could not fix fn_datafolder, gonna skip"
						print "Final fn_datafolder", emb.fn_datafolder
						raw_input()
						continue
					
					if not os.path.isfile(emb.fn_preimage):
						print "Warning, really could not fix fn_preimage, gonna skip"
						print "Final fn_preimage", emb.fn_preimage
						raw_input()
						continue
				
			#Reload file list for good measure
			emb.update_file_list()
			
			#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#Select boundary
			
			if bound:
				selector=pyfrp_plot.boundary_selector(emb)
				emb=selector.get_embryo()
			
			
			#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#Analyze embryo

			if analyze==1:
				
				emb = pyfrp_img.analyze_dataset(emb)	
				print "Done Image Analysis"
				#emb.plot_data()
				mol.save_molecule(fn_out_mol)	
				
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
				
				#emb.conc_rim=emb.out_av_data_d[0]
				
				emb.conv_fact=566.79/emb.data_res_px
				
				#Some settings
				#Plots
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
				emb.volSize_px=volSize_px
				
				if geometry==0:
					emb.geometry="Cylinder"
					
					if emb.quad_red:
						emb.fn_mesh="meshfiles/quad_cyl"
						fn_geo="meshfiles/quad_cyl"
					else:
						emb.fn_mesh="meshfiles/cylinder"
						fn_geo="meshfiles/cylinder"
					
					
					if refine_slice:

						emb.fn_mesh=emb.fn_mesh+"_field_box"
						fn_geo=fn_geo+"_field_box"
					
					emb.fn_mesh=emb.fn_mesh+".msh"
					fn_geo=fn_geo+".geo"
					
				
					
					if refine_slice:
						#update_refined_cylinder_msh(fn,radius,height,center,volSize_px,sl_height,sl_width,sl_extend,vol,run=True,debug=False)					
						
						pyfrp_gmsh.update_refined_cylinder_msh(fn_geo,emb.radius_embr_px,emb.cylinder_height_px,emb.center_embr_px,emb.volSize_px,emb.slice_depth_px,emb.slice_width_px,emb.radius_embr_px,emb.volSize_px/2.)
					else:
						pyfrp_gmsh.update_cylinder_msh(fn_geo,emb.radius_embr_px,emb.cylinder_height_px,emb.center_embr_px,emb.volSize_px)
					
					
					pyfrp_gmsh.refine_msh(emb.fn_mesh)
				
				elif geometry==1:
					emb.geometry="Fish"
					
					if emb.quad_red:
						emb.fn_mesh="meshfiles/quad_dome"
						fn_geo="meshfiles/quad_dome"
					else:
						emb.fn_mesh="meshfiles/dome"
						fn_geo="meshfiles/dome"
					
					pyfrp_gmsh.update_dome_msh(fn_geo,emb.radius_embr_px,-emb.slice_depth_px,emb.center_embr_px,emb.volSize_px)
					pyfrp_gmsh.refine_msh(emb.fn_mesh)
					
				#Make sure to average everything for control
				emb.avg_all=1
				
				#Debug?
				emb.debug_simulation=1
					
				#Apply ICs via interpolation
				emb.apply_data=3
				
				#Reference diffusion coefficient
				emb.D=10.
				
				#Timesteps
				emb.steps_sim=steps_sim
					
				#Converty to logarithmic scale
				#emb.to_log_scale()
				emb.to_lin_scale()
				emb=pyfrp_fit.get_opt_tvec_sim_px(emb,150)
				
				#run PDE solver
				print emb.fn_mesh
				print fn_geo
				emb=pyfrp_sim.simulate_diff_react(emb)
				
				mol.save_molecule(fn_out_mol)
			
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
				emb.add_fit(0,"")
				
				emb.fits[0].opt_meth="Constrained Nelder-Mead"
				emb.fits[0].fit_prod=0
				emb.fits[0].fit_degr=0
				emb.fits[0].fit_pinned=1
				emb.fits[0].fit_slice=fit_slice
				emb.fits[0].fit_squ=fit_squ
				emb.fits[0].fit_out=fit_out
				emb.fits[0].UB_D=150.
				emb.fits[0].LB_D=0.01
				emb.fits[0].x0=[10,0,0]
				emb.fits[0].equ_on=1
				emb.fits[0].fit_cut_off_t=0
				
				#Fit unpinned
				emb.add_fit(1,"")			

				emb.fits[1].opt_meth="Constrained Nelder-Mead"
				emb.fits[1].fit_prod=0
				emb.fits[1].fit_degr=0
				emb.fits[1].fit_slice=fit_slice
				emb.fits[1].fit_squ=fit_squ
				emb.fits[1].fit_out=fit_out
				emb.fits[1].fit_pinned=0
				emb.fits[1].UB_D=150.
				emb.fits[0].LB_D=0.01
				emb.fits[1].x0=[10,0,0]
				emb.fits[1].equ_on=0
				
				#Fit equ, pinned, cut_off
				emb.add_fit(2,"")			

				emb.fits[2].opt_meth="Constrained Nelder-Mead"
				emb.fits[2].fit_prod=0
				emb.fits[2].fit_degr=0
				emb.fits[2].fit_pinned=1
				emb.fits[2].fit_slice=fit_slice
				emb.fits[2].fit_squ=fit_squ
				emb.fits[2].fit_out=fit_out
				emb.fits[2].UB_D=400
				emb.fits[0].LB_D=0.01
				emb.fits[2].x0=[10,0,0]
				emb.fits[2].equ_on=1
				emb.fits[2].fit_cut_off_t=1
				emb.fits[2].cut_off_t=50
				
				#Fit cut_off
				emb.add_fit(3,"")			

				emb.fits[3].opt_meth="Constrained Nelder-Mead"
				emb.fits[3].fit_prod=0
				emb.fits[3].fit_degr=0
				emb.fits[3].fit_slice=fit_slice
				emb.fits[3].fit_squ=fit_squ
				emb.fits[3].fit_out=fit_out
				emb.fits[3].fit_pinned=0
				emb.fits[3].UB_D=400
				emb.fits[0].LB_D=0.01
				emb.fits[3].x0=[10,0,0]
				emb.fits[3].equ_on=0
				emb.fits[3].fit_cut_off_t=1
				emb.fits[3].cut_off_t=50
				
				#Fit equ, pinned TNC
				emb.add_fit(4,"")
				
				emb.fits[4].opt_meth="TNC"
				emb.fits[4].fit_prod=0
				emb.fits[4].fit_degr=0
				emb.fits[4].fit_pinned=1
				emb.fits[4].fit_slice=fit_slice
				emb.fits[4].fit_squ=fit_squ
				emb.fits[4].fit_out=fit_out
				emb.fits[4].UB_D=95.
				emb.fits[4].LB_D=0.01
				emb.fits[4].x0=[10,0,0]
				emb.fits[4].equ_on=1
				emb.fits[4].fit_cut_off_t=0
				
				#Perform fits
				emb=pyfrp_fit.parm_fitting(emb,0)
				emb=pyfrp_fit.parm_fitting(emb,1)
				#emb=pyfrp_fit.parm_fitting(emb,2)
				#emb=pyfrp_fit.parm_fitting(emb,3)
				emb=pyfrp_fit.parm_fitting(emb,4)
				
				#Add pinned fit to selected fit
				mol.sel_fits.append(emb.fits[0])
				
				print "Done Fitting"
				mol.save_molecule(fn_out_mol)
				
				
			if out:
				
				print "-----------------------"
				print "pinned eqon CNM:"
				emb.fits[0].print_results()
				emb.fits[0].plot_fit_pinned()
				print "-----------------------"
				print "unpinned eqoff:"
				emb.fits[1].print_results()
				emb.fits[1].plot_fit_unpinned()
				#print "-----------------------"
				#print "pinned eqon cut_off:"
				#emb.fits[2].print_results()
				#emb.fits[2].plot_fit_pinned()
				#print "unpinned eqoff cut_off:"
				#emb.fits[3].print_results()
				#emb.fits[3].plot_fit_unpinned()
				#print "-----------------------"
				#print "pinned eqon TNC:"
				#emb.fits[4].print_results()
				#emb.fits[4].plot_fit_pinned()
				raw_input()
	print "Saved molecule to: ", fn_out_mol
	
	
	if stats:
		
		print "-------------------------------------------------------------------"
		print "Final statistics:"
		print "-------------------------------------------------------------------"
		
		
		mol.sumup_results()
		
		print "mol.D_mu_av=",mol.D_mu_av
		print "mol.degr_av=",mol.degr_av
		print "mol.prod_av=",mol.prod_av
		print "mol.Rsq_av=",mol.Rsq_av
		print "mol.ssd_av=",mol.ssd_av
		
		print "mol.D_mu_std=",mol.D_mu_std
		print "mol.prod_std=",mol.prod_std	
		print "mol.degr_std=",mol.degr_std
		
		print "mol.D_mu_sterr=",mol.D_mu_sterr
		print "mol.prod_sterr=",mol.prod_sterr		
		print "mol.degr_sterr=",mol.degr_sterr
		
		raw_input()
