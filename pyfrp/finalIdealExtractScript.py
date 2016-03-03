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
print "Will extract embryos of molecule files in the folder  " + fn + " to " + fnOut 

print "================================================================================"

#===========================================================================================================================================================================
#Script
#===========================================================================================================================================================================

#Get all molecule files
molFiles=pyfrp_misc_module.getSortedFileList(fn,".mol")

for f in molFiles:
	
	mol=pyfrp_IO_module.loadFromPickle(fn+f)
	
	try:
		os.system("mkdir -v "+ fn+fnOut)
	except OSError:
		pass
	
	try:
		os.system("mkdir -v " + fn+fnOut + mol.name)
	except OSError:
		pass
	
	mol.extractEmbryos2Files(fn=fn+fnOut+mol.name+"/")
	
	print "Done with molecule "+ f 
	
print "Done"	
	