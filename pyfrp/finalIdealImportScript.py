###Final script to extract embryo files from molecule files for IdealFRAP Project

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#PyFRAP
import pyfrp_IO_module
import pyfrp_misc_module

#Misc
import os
from optparse import OptionParser

#===========================================================================================================================================================================
#Defining parsing possibilities
#===========================================================================================================================================================================

usage = "usage: %prog [options]"
parser = OptionParser(usage)
     
parser.add_option("-f", "--fn",action="store", type='string', dest="fn", help="Folder containing molecule data (eg.FITC-dextran-10kDa)")
parser.add_option("-o", "--out",action="store", type='string', dest="fn", help="Folder containing ebmryo files data (eg.embryoFiles/)")

#Set defaults 
parser.set_defaults(fn=None)
parser.set_defaults(out="embryoFiles/")

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
fnOut=str(options.out)

print "================================================================================"
print "Will import all embryos for the molecule files in the folder  " + fn + " to " + fnOut 

print "================================================================================"

#===========================================================================================================================================================================
#Script
#===========================================================================================================================================================================

#Get all molecule files
molFiles=pyfrp_misc_module.getSortedFileList(fn,".mol")

#Loop through all molecule files
for f in molFiles:
	
	#Load molecule
	mol=pyfrp_IO_module.loadFromPickle(fn+f)
	
	#Get all embryo files
	embFiles=pyfrp_misc_module.getSortedFileList(fn+fnOut+mol.name,".emb")
	
	#Loop through embryos
	for i,emb in enumerate(mol.embryos):
		
		#Loop through all embryo files
		for embFn in embFiles:
			embFnStripped=embFn.replace(".emb","")
			if embFnStripped==emb.name:
				
				#Remove embryo
				mol.embryos.pop(i)
				
				#Load new embryo
				newEmbryo=pyfrp_IO_module.loadFromPickle(fn+fnOut+mol.name+"/"+embFn)
				
				#Insert it in embryos list
				mol.embryos.insert(i,newEmbryo)
				
				print "Successfully imported ", embFn
				
				break
			
	#Save molecule again
	mol.save(fn+f)
	
	print "Done with molecule "+ f 
	
	
	
print "Done"	
	