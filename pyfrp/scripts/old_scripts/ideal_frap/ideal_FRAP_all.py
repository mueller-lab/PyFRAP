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
from molecule import *
import time


#Some flags what to do
analyze=0
preproc=0
simulate=0
pin=0
fit=0
show_results=1
build_mol=0
sum_datasets=0
save_plots=0
gauss=0

#Define lists of datasets

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#June A488 Datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#datasets=["A488-dextran-3kDa_20140612","A488-dextran-3kDa_1uM_agarose_1Percent_20140616","A488-dextran-10kDa_20140610","A488-dextran-10kDa_agarose_1percent_20140616"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#June/July FITC Datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#datasets=["fluorescein-dextran-3kDa_1uM_H2O_20140701"]
#datasets=["fluorescein-dextran-4kDa_sigma_100nM_H2O_20140710"]
#datasets=["FITC-dextran-70kDa_sigma_1uM_H2O_20140715"]
#datasets=["FITC-dextran-70kDa_sigma_1uM_H2O_20140714"]
#datasets=["fluorescein-dextran-40kDa_1uM_H2_20140630"]
#datasets=["FITC-dextran-150kDa_sigma_1uM_H2O_20140715"]
#exps_all=[["exp1"]]
#exps_all=[["exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10"]]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#August FITC Datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#datasets=["FITC-dextran-150kDa_sigma_1uM_H2O_20140805","FITC-dextran-70kDa_sigma_1uM_H2O_20140807","FITC-dextran-70kDa_sigma_1uM_H2O_20140807"]
#exps_all=[["exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10","exp11","exp12"],["exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10"]]

#datasets=["FITC-dextran-500kDa_1uM_H2O_20140821","FITC-dextran-10kDa_10uM_H2O_20140821","FITC-dextran-4kDa_sigma_1uM_H2O_20140801","FITC-dextran-3kDa_invitrogen_1uM_H2O_20140805","FITC-dextran-40kDa_invitrogen_1uM_H2O_20140801"]
#exps_all=[["exp1","exp3","exp4","exp5","exp6"],["exp2","exp3"],["exp9","exp10","exp11","exp12"],["exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8","exp9"],["exp1","exp2","exp3","exp4","exp5","exp6","exp7","exp8"]]

#datasets=["FITC-dextran-500kDa_1uM_H2O_20140828","FITC-dextran-10kDa_10uM_H2O_20140828"]#
#exps_all=[["exp1","exp2","exp3"],["exp1","exp2","exp4","exp5","exp6","exp7","exp8","exp9","exp10"]]

datasets=["FITC-dextran-70kDa_invitrogen_1uM_H2O_20140825","FITC-dextran-70kDa_sigma_1uM_H2O_20140731"]
exps_all=[["exp1","exp3","exp4","exp5","exp6","exp7","exp8","exp9","exp10"],["exp1","exp2","exp3","exp4","exp5"]]

#"exp6","exp7","exp8","exp9","exp10"]]

#datasets=["FITC-dextran-70kDa_sigma_1uM_H2O_20140731"]
#exps_all=[["exp7","exp8","exp9","exp10"]]


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating embryo object

j=0
#oldDs=[]
#newDs=[]

for dataset in datasets:
	
	exps=exps_all[j]
	j=j+1
	
	if build_mol==1:
		mol=molecule(dataset)
	
	for exp in exps:
		
		print "======================================================================================================================"
		print "Now processing", dataset, exp
		print "======================================================================================================================"
		
		#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#Load embryo object
		
		emb=embryo("blub","frap")
		emb.name=dataset+"_"+exp
		
		#Updating some folders
		emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
		
		
		fn_load=emb.fn_resultfolder+"/"+emb.name+".pk"
		emb=emb.load_embryo(fn_load)	
		
		emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
		emb.fn_datafolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/post"
		
		prefolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/pre/"
		prelist=get_sorted_folder_list(prefolder,emb.data_ft)
		emb.fn_preimage=prefolder+prelist[0]
		
		fn_save=emb.fn_resultfolder+"/"+emb.name+"_gauss.pk"
		
		emb.fn_mesh=emb.fn_resultfolder+"/meshfile.msh"
		
		
		
		#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#Analyze embryo

		if analyze==1:
			
			
			emb.gaussian=gauss
			emb.gaussian_sigma=3
			
			emb.file_list=get_sorted_folder_list(emb.fn_datafolder,emb.data_ft)
			
			print emb.file_list[0]
			emb = analyze_dataset(emb)
			if gauss==1:
				emb.save_embryo(fn_save)
			else:	
				emb.save_embryo(None)
			
			print "Done Image Analysis"
	
		#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#Preproc for ICs 
			
		if preproc==1:
			
			emb = apply_data_ics(emb)
			if gauss==1:
				emb.save_embryo(fn_save)
			else:	
				emb.save_embryo(None)
			
			
			print "Done Preproc"
			
		#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		#Simulate for scaling solution 

		if simulate==1:
			
			#Load mesh
			emb.mesh=GmshImporter3D(emb.fn_mesh)
			emb.steps_sim=3000
			emb=simulate_diff_react(emb)
			if gauss==1:
				emb.save_embryo(fn_save)
			else:	
				emb.save_embryo(None)
			
			print "Done Simualtion"

			
		if pin==1:
			emb.debug_pinning=0
			emb=pin_conc(emb)	
			if gauss==1:
				emb.save_embryo(fn_save)
			else:	
				emb.save_embryo(None)
			
			print "Done Pinning"
			
			
		if fit==1:
			
			#Deleting all fits
			emb.fits=[]
			emb.debug_fit=0
			
			#Fit equ, pinned
			emb.add_fit(0,"fit equ pinned","default")

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
			emb.add_fit(1,"fit unpinned","default")			

			emb.fits[1].opt_meth="Nelder-Mead"
			emb.fits[1].fit_prod=0
			emb.fits[1].fit_degr=0
			emb.fits[1].fit_slice=0
			emb.fits[1].fit_squ=1
			emb.fits[1].fit_out=0
			emb.fits[1].fit_pinned=0
			emb.fits[1].UB_D=400
			emb.fits[1].x0=[50,0,0]
			emb.fits[1].equ_on=0
		
		
			#Fit equ, pinned, cut_off
			emb.add_fit(2,"fit equ pinned cut_off","default")			

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
			emb.add_fit(3,"fit cut_off","default")			

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
			
			#Fit equ, pinned, cut_off long
			emb.add_fit(4,"fit equ pinned cut_off long","default")

			emb.fits[4].opt_meth="Nelder-Mead"
			emb.fits[4].fit_prod=0
			emb.fits[4].fit_degr=0
			emb.fits[4].fit_pinned=1
			emb.fits[4].fit_slice=0
			emb.fits[4].fit_squ=1
			emb.fits[4].fit_out=0
			emb.fits[4].UB_D=400
			emb.fits[4].x0=[50,0,0]
			emb.fits[4].equ_on=1
			emb.fits[4].fit_cut_off_t=1
			emb.fits[4].cut_off_t=200
			
			
			#Fit unpinned, cut_off long
			emb.add_fit(5,"Fit unpinned cut_off long","default")			

			emb.fits[5].opt_meth="Nelder-Mead"
			emb.fits[5].fit_prod=0
			emb.fits[5].fit_degr=0
			emb.fits[5].fit_slice=0
			emb.fits[5].fit_squ=1
			emb.fits[5].fit_out=0
			emb.fits[5].fit_pinned=0
			emb.fits[5].UB_D=400
			emb.fits[5].x0=[50,0,0]
			emb.fits[5].equ_on=0
			emb.fits[5].fit_cut_off_t=1
			emb.fits[5].cut_off_t=200
			
			#Perform fits
			emb=parm_fitting(emb,0)
			emb=parm_fitting(emb,1)
			emb=parm_fitting(emb,2)
			emb=parm_fitting(emb,3)
			emb=parm_fitting(emb,4)
			emb=parm_fitting(emb,5)
			
			if gauss==1:
				emb.save_embryo(fn_save)
			else:	
				emb.save_embryo(None)
			

			print "Done Fitting"
			#newDs.append(emb.fits[0].D_opt_mu)
			
			fn_save=emb.fn_resultfolder+"/"+"resultcsv.csv"
			
			wfile=csv.writer(open(fn_save,'wb'), delimiter=';')
			
			names=["fit_pinned","equ_on","fit_cut_off_t","cut_off_t","fit_squ","fit_slice","fit_out","D_opt_mu (um^2/s)","success","ssd","Rsq","Comments"]
			
			wfile.writerow(names)
			
			for thefit in emb.fits:
			
				vals=[thefit.fit_pinned,thefit.equ_on,thefit.fit_cut_off_t,thefit.cut_off_t,thefit.fit_squ,thefit.fit_slice,thefit.fit_out,thefit.D_opt_mu,thefit.success,thefit.ssd,thefit.Rsq,""]
			
				wfile.writerow(vals)
				
		
		if show_results==1:
			
			
			emb.name=exp
			
			print "-----------------------"
			print "pinned eqon:"
			emb.fits[0].print_results()
			emb.fits[0].plot_fit_pinned()
			
			
			#print "-----------------------"
			#print "unpinned eqoff:"
			#emb.fits[1].print_results()
			#emb.fits[1].plot_fit_unpinned()
			
			
			#print "-----------------------"
			#print "pinned eqon cut_off:"
			#emb.fits[2].print_results()
			#emb.fits[2].plot_fit_pinned()
			#print "-----------------------"
			#print "unpinned eqoff cut_off:"
			#emb.fits[3].print_results()
			#emb.fits[3].plot_fit_unpinned()
			#print "-----------------------"
			#print "pinned eqon cut_off long:"
			#emb.fits[4].print_results()
			#emb.fits[4].plot_fit_pinned()
			#print "-----------------------"
			#print "unpinned eqoff cut_off long:"
			#emb.fits[5].print_results()
			#emb.fits[5].plot_fit_unpinned()
			
			
		if build_mol==1:
			mol.add_embryo(emb)
		
		if save_plots==1:
			for fit in emb.fits:
				fn_save="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/results/figs/"+emb.name+"_"+fit.name+".png"
				
				if fit.fit_pinned==1:
					fit.save_plot_fit_pinned(fn_save)
				else: 
					fit.save_plot_fit_unpinned(fn_save)
				
			
		
	if build_mol==1:
		
		#Write molecule results into csv file
		fn_save="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/results/"+mol.name+".csv"
			
		wfile=csv.writer(open(fn_save,'wb'), delimiter=';')
		
		names=["dataset","exp","fit_pinned","equ_on","fit_cut_off_t","cut_off_t","fit_squ","fit_slice","fit_out","D_opt_mu (um^2/s)","success","ssd","Rsq","Comments"]
			
		wfile.writerow(names)
		
		fit_count=len(mol.embryos[0].fits)
		
		for fit_num in range(fit_count):
			
			mol.sel_fits=[]
			
			for emb in mol.embryos:
			
				mol.sel_fits.append(emb.fits[fit_num])
				
				
			mol.sumup_results()
			
			print "-----------------------"
			print "Molecule: ", mol.name
			print "Type:", mol.sel_fits[0].name
			print "Average D_mu_av=:", mol.D_mu_av
			print "D_mu_std=:", mol.D_mu_std
			print "D_mu_sterr=:", mol.D_mu_sterr
			
			
			wfile.writerow([mol.sel_fits[0].name])
			
			for thefit in mol.sel_fits:
			
				vals=[dataset,thefit.embryo.name,thefit.fit_pinned,thefit.equ_on,thefit.fit_cut_off_t,thefit.cut_off_t,thefit.fit_squ,thefit.fit_slice,thefit.fit_out,thefit.D_opt_mu,thefit.success,thefit.ssd,thefit.Rsq,""]
				
				wfile.writerow(vals)
			
			wfile.writerow([])
			wfile.writerow(["Averaged:"])
			wfile.writerow([dataset,"",thefit.fit_pinned,thefit.equ_on,thefit.fit_cut_off_t,thefit.cut_off_t,thefit.fit_squ,thefit.fit_slice,thefit.fit_out,mol.D_mu_av,"",mol.ssd_av,mol.Rsq_av,""])
			wfile.writerow([])
			wfile.writerow([])
			
		mol.save_molecule("/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/results/"+dataset+".pk")
		
if sum_datasets==1:
	
	#Write molecule results into csv file
	fn_save="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/summary_june.csv"
			
	wfile=csv.writer(open(fn_save,'wb'), delimiter=';')
		
	names=["dataset","","fit_pinned","equ_on","fit_cut_off_t","cut_off_t","fit_squ","fit_slice","fit_out","D_opt_mu (um^2/s)","","ssd","Rsq","Comments"]
	
	
	for fit_num in range(fit_count):
		
		count=0
		for dataset in datasets:
	
			mol=molecule("blub")
			mol=mol.load_molecule("/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/results/"+dataset+".pk")
				
			mol.sel_fits=[]
			
			for emb in mol.embryos:
				mol.sel_fits.append(emb.fits[fit_num])
				
			mol.sumup_results()
			
			if count==0:
				wfile.writerow([mol.sel_fits[0].name])
				print "Writing summary of all fits with property:", mol.sel_fits[0].name 
			print "Writing dataset:", dataset 
				
			thefit=mol.sel_fits[0]
			wfile.writerow([dataset,"",thefit.fit_pinned,thefit.equ_on,thefit.fit_cut_off_t,thefit.cut_off_t,thefit.fit_squ,thefit.fit_slice,thefit.fit_out,mol.D_mu_av,"",mol.ssd_av,mol.Rsq_av,""])	
			
			count=count+1
			
		wfile.writerow([])
		wfile.writerow([])

#for i in range(len(oldDs)):
	#print oldDs[i], newDs[i], abs(oldDs[i]-newDs[i])/oldDs[i]
	
