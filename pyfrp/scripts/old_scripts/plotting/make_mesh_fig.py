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



#Choosing dataset

emb=embryo("blub","default")


fn_load="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/myembryo.pk"
emb=emb.load_embryo(fn_load)
emb.fn_resultfolder="/home/alex_loc/Documents/Research/PyFRAP/Test_datasets/"+"011912_Lft1-LGDPPVAT-GFP"+"/"+"emb4"+"/results/"
emb.fn_mesh=emb.fn_resultfolder+"meshfile_cyl.msh"
emb.mesh=GmshImporter3D(emb.fn_mesh)

#USAGE: plot3Dmesh(mesh,region,high_region,bound_region,alpha,markers)

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#region options
#1 -> plot whole mesh (blue)
#2 -> plot only boundary vertices, can be restricted to bound_region, can be either given as a list of node indices or [[xmin,xmax],yxmin,ymax],[zmin,zmax]](green)
#3 -> plot only highlighted region, high_region can be either given as a list of node indices or as a list [[xmin,xmax],yxmin,ymax],[zmin,zmax]] (red)
#4 -> plot highlighted+whole mesh, mesh transparent
#5 -> plot highlighted+boundary, boundary transparent
#alpha = degree of transparency, 0=totally transparent, 1=not transparent
#----------------------------------------------------------------------------------------------------------------------------------------------------------


#Crosssection
#plot3Dmesh(emb.mesh,2,[[256-100,256+100],[256-100,256+100],[":",":"]],[[256,":"],[":",":"],[":",":"]],0.3,'')

#Just boundary
plot3Dmesh(emb.mesh,2,[],[],0.3,'')

#Boundary + slice

#plot3Dmesh(emb.mesh,5,[[":",":"],[":",":"],[-60,-50]],[],0.3,'')

