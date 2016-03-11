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
     
parser.add_option("-n", "--name",action="store", type='string', dest="n", help="Name of experiment")
parser.add_option("-d", "--date",action="store", type='string', dest="d", help="Date of experiment")
parser.add_option("-f", "--folder",action="store", type='string', dest="f", help="Data folder")
parser.add_option("-B", "--beads",action="store", type='string', dest="B", help="Beads molecule file")
parser.add_option("-N", "--nobeads",action="store", type='string', dest="N", help="NoBeads molecule file")
parser.add_option("-E", "--embryo",action="store", type='string', dest="E", help="Embryo molecule file")

parser.set_defaults(d="")
parser.set_defaults(f="")
parser.set_defaults(n="")
parser.set_defaults(B="")
parser.set_defaults(N="")
parser.set_defaults(E="")

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.f)
date=str(options.d)
name=str(options.n)
fn_beads=str(options.B)
fn_nobeads=str(options.N)
fn_embryo=str(options.E)

print "================================================================================"
print "Will sort the following folder" + fn
print "name: ", name
print "date: ", date 
print "beads: ", fn_beads 
print "nobeads: ", fn_nobeads
print "embryo: ", fn_embryo
print "================================================================================"

#===========================================================================================================================================================================
#Preperations
#===========================================================================================================================================================================

beads_exists=False
nobeads_exists=False
embryo_exists=False

if os.path.isfile(fn+"/"+fn_beads):
	beads_exists=True
	mol_beads=molecule("blub")
	mol_beads=mol_beads.load_molecule(fn+"/"+fn_beads)
else:
	print "Warning, " +fn_beads+  " does not exist"

if os.path.isfile(fn+"/"+fn_nobeads):
	nobeads_exists=True
	mol_nobeads=molecule("blub")
	mol_nobeads=mol_nobeads.load_molecule(fn+"/"+fn_nobeads)
else:
	print "Warning, " +fn_nobeads+  " does not exist"
	
if os.path.isfile(fn+"/"+fn_embryo):
	embryo_exists=True
	mol_embryo=molecule("blub")
	mol_embryo=mol_embryo.load_molecule(fn+"/"+fn_embryo)
else:
	print "Warning, " +fn_embryo+  " does not exist"	
	
#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 0: Do some checks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 0: Do some checks
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#Grab path of folder
root,fn_name=os.path.split(fn)

#Is folder with name=name already there? Then use that one
print "Checking if ", root+name, "already exists:" 
if os.path.isdir(root+"/"+name):
	print "It already exists. Going to change fn to it!"
	fn_dest=root+"/"+name
else:
	print "It does not already exists."
	fn_dest=fn
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Make folder structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Make folder structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
try:
	os.mkdir(fn_dest+"/beads")
except OSError:	
	print "Could not make folders:"
	print fn_dest+"/beads"
	
try:
	os.mkdir(fn_dest+"/nobeads")
except OSError:	
	print "Could not make folders:"
	print fn_dest+"/nobeads"
	
try:	
	os.mkdir(fn_dest+"/embryo")
except OSError:	
	print "Could not make folders:"
	print fn_dest+"/embryo"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Move and rename embryo folders name
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 2: Move and rename embryo folders name
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

files=os.listdir(fn)
for f in files:
	if os.path.isdir(fn+"/"+f):
		if f not in ["lsm","beads","nobeads"]:
			print f
			for emb in mol_embryo.embryos:
				print emb.name, emb.fn_datafolder
				print 
			
			if f in ["beads","nobeads","embryo"]:
				pass
			else:		
				if "_nobeads" in f:
					
					#Grab embryo
					emb=mol_nobeads.emb_by_name("_"+str(int(f.strip("_nobeads"))))
					
					if emb!=False:
						
						emb.fn_datafolder=pyfrp_misc.win2lin_path(emb.fn_datafolder)
						
						#Update fn_datafolder
						
						if emb.fn_datafolder.endswith("/"):
							head,tail=os.path.split(emb.fn_datafolder[:-1])
						else:
							head,tail=os.path.split(emb.fn_datafolder)
						emb.fn_datafolder=root+"/"+name+"/nobeads/"+date+"_"+f.strip("_nobeads")+"/"+tail+"/"
					
						#Update fn_preimage
						head,tail=os.path.split(emb.fn_preimage)
						emb.fn_preimage=root+"/"+name+"/nobeads/"+date+"_"+f.strip("_nobeads")+"/"+"pre/"+tail
						
						#Move and rename
						new_fn=fn_dest+"/nobeads/"+date+"_"+f.strip("_nobeads")
						
						os.system("mv -v " + fn+"/"+f+" "+new_fn)
						
				elif "_embryo" in f:
					
					#Grab embryo
					emb=mol_embryo.emb_by_name("_"+str(int(f.strip("_embryo"))))
					
					if emb!=False:
						
						emb.fn_datafolder=pyfrp_misc.win2lin_path(emb.fn_datafolder)
						
						#Update fn_datafolder
						if emb.fn_datafolder.endswith("/"):
							head,tail=os.path.split(emb.fn_datafolder[:-1])
						else:
							head,tail=os.path.split(emb.fn_datafolder)
						emb.fn_datafolder=root+"/"+name+"/embryo/"+date+"_"+f.strip("_embryo")+"/"+tail+"/"
					
						#Update fn_preimage
						head,tail=os.path.split(emb.fn_preimage)
						emb.fn_preimage=root+"/"+name+"/embryo/"+date+"_"+f.strip("_embryo")+"/"+"pre/"+tail
						
						#Move and rename
						new_fn=fn_dest+"/embryo/"+date+"_"+f.strip("_embryo")
						
						os.system("mv -v " + fn+"/"+f+" "+new_fn)
			
				else:
					
					#Grab embryo
					emb=mol_beads.emb_by_name("_"+str(int(f.strip("_nobeads"))))
					
					if emb!=False:
						#Update fn_datafolder
						if emb.fn_datafolder.endswith("/"):
							head,tail=os.path.split(emb.fn_datafolder[:-1])
						else:
							head,tail=os.path.split(emb.fn_datafolder)
						emb.fn_datafolder=root+"/"+name+"/beads/"+date+"_"+f.strip("_nobeads")+"/"+tail
						
						#Update fn_preimage
						head,tail=os.path.split(emb.fn_preimage)
						emb.fn_preimage=root+"/"+name+"/beads/"+date+"_"+f.strip("_nobeads")+"/"+"pre/"+tail
						
						#Move and rename
						new_fn=fn_dest+"/beads/"+date+"_"+f
						
						os.system("mv -v " + fn+"/"+f+" "+new_fn)
				
				#Rename embryo
				if emb!=False:
					emb.name=date+"_emb"+f
				else:
					print "Warning, could not find embryo containing:" + "_"+str(int(f.strip("_nobeads"))) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Set name molecule file and save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 3: Set name of molecule file and save
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

if beads_exists:
	mol_beads.name=name+"_beads_"+date
	fn_beads_new=fn+"/"+name+"_beads_"+date+".pk"
	mol_beads.save_molecule(fn_beads_new)

if nobeads_exists:
	mol_nobeads.name=name+"_nobeads_"+date
	fn_nobeads_new=fn+"/"+name+"_nobeads_"+date+".pk"
	mol_nobeads.save_molecule(fn_nobeads_new)

if embryo_exists:
	mol_embryo.name=name+"_embryo_"+date
	fn_embryo_new=fn+"/"+name+"_embryo_"+date+".pk"
	mol_embryo.save_molecule(fn_embryo_new)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 4: Rename molecule file and move it if necessary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 4: Rename molecule file and move it if necessary
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

if beads_exists:
	os.system("mv -v "+ fn_beads_new+" "+ fn_dest+"/"+name+"_beads_"+date+".pk")
if nobeads_exists:
	os.system("mv -v "+ fn_nobeads_new+" "+ fn_dest+"/"+name+"_nobeads_"+date+".pk")
if embryo_exists:
	os.system("mv -v "+ fn_embryo_new+" "+ fn_dest+"/"+name+"_embryo_"+date+".pk")
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 5: Rename main folder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 5: Rename main folder
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

print "mv -v "+fn_dest+" " + root+"/"+name
os.system("mv -v "+fn_dest+" " + root+"/"+name)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 6: Clean up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 6: Clean Up
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

files=os.listdir(root+"/"+name)
for f in files:
	if f.endswith(".pk") and "Analysis" in f:
		os.system("rm -v " + root+"/"+name+"/"+f)
