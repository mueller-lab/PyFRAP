#===========================================================================================================================================================================
#Importing necessary modules
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
parser.add_option("-O", "--fnout",action="store", type='string', dest="O", help="Filename for output")
parser.add_option("-n", "--name",action="store", type='string', dest="n", help="Name for output")


parser.set_defaults(fn=None)

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
fn_out=str(options.O)
name=str(options.n)


print "================================================================================"
print "Will group all datasets specified in file " + fn
print "Will do the following:"
print "fn_out: ", fn_out
print "name: ", name
print "================================================================================"

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
#Run script
#===========================================================================================================================================================================

#Create master molecule file
mol_all=molecule(name)

#Loop through all molecules specified
for k,fn_mol in enumerate(mols):
	
	print "======================================================================================================================================="
	print "Processing molecule " + fn_mol
	print "======================================================================================================================================="
	
	#Creating molecule object and load molecule file
	mol=molecule("blub")
	mol=mol.load_molecule(fn_mol)
	
	#Loop through embryos and append to master molecule file
	for emb in mol.embryos:
		mol_all.embryos.append(emb)
	
	print "Appended ", len(mol.embryos), " embryos to molecule ", name 
	
#Save master molecule file	
mol_all.save_molecule(fn_out)