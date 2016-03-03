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
from pyfrp_gmsh_module import *

from embryo import *
import time

#Define dataset
dataset="FITC-dextran-70kDa_invitrogen_1uM_H2O_20140825"

exps_num=arange(1,11)

exps_all=[]
for i in exps_num:

	exps_all.append("exp"+str(int(i)))



for exp in exps_all:

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Creating embryo object
	emb=embryo("blub","frap")
	emb.name=dataset+"_"+exp
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
	
	fn_load=emb.fn_resultfolder+"/"+emb.name+".pk"
	emb=emb.load_embryo(fn_load)
	
	#Dataset
	emb.fn_datafolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/post"
	emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/results"
		
	if emb.norm_by_pre==1:
		prefolder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/"+exp+"/pre/"
		prelist=get_sorted_folder_list(prefolder,emb.data_ft)
		emb.fn_preimage=prefolder+prelist[0]
		
	emb.fn_mesh=emb.fn_resultfolder+"/meshfile.msh"
	emb.fn_geo=emb.fn_resultfolder+"/meshfile.geo"
	
	emb.save_embryo(None)
	print "Changed filepaths of embryo: ", emb.name
	
	
	