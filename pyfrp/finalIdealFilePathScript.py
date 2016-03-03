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

#Set defaults 
parser.set_defaults(fn=None)

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)


print "================================================================================"
print "Will change filepaths of embryos of molecule files in the folder  " + fn

print "================================================================================"

#===========================================================================================================================================================================
#Script
#===========================================================================================================================================================================

#Define paths
basepath='../data/Gary'
splitpath='Gary'

#Get all molecule files
molFiles=pyfrp_misc_module.getSortedFileList(fn,".mol")

for f in molFiles:
	
	#Load molecule
	mol=pyfrp_IO_module.loadFromPickle(fn+f)
	
	#Loop through all embryos and reset filepaths
	for emb in mol.embryos:
		datafolder=emb.getDataFolder()
		dump,datafolder=datafolder.split(splitpath)
		datafolder=basepath+datafolder
		
		emb.setDataFolder(datafolder)
		
		if 'post' in datafolder:
			prefolder=datafolder.replace('post','pre')
		else:	
			prefolder=datafolder.replace('recover','pre')
		
		emb.analysis.setFnPre(prefolder)
		emb.analysis.setFnFlatten('../data/flatteningData/20152011_flat/01/')
		emb.analysis.setFnBkgd('../data/backgroundData/20160105_bkgd/01/')
		
	#Save molecule file again
	mol.save(fn+f)
		
	print "Done with molecule "+ f 
	
print "Done"	
	