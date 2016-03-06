###Final script to extract embryo files from molecule files for IdealFRAP Project

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#PyFRAP
import pyfrp_IO_module
import pyfrp_misc_module
import pyfrp_plot_module

#Misc
import os
from optparse import OptionParser

#Matplotlib
import matplotlib.pyplot as plt

#numpy
import numpy as np

#===========================================================================================================================================================================
#Defining parsing possibilities
#===========================================================================================================================================================================

usage = "usage: %prog [options]"
parser = OptionParser(usage)
     
parser.add_option("-f", "--fn",action="store", type='string', dest="fn", help="Folder containing molecule data (eg.FITC-dextran-10kDa)")
parser.add_option("-o", "--out",action="store", type='string', dest="fn", help="Folder containing ebmryo files data (eg.results/)")

#Set defaults 
parser.set_defaults(fn=None)
parser.set_defaults(out="results/")

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
fnOut=str(options.out)

print "================================================================================"
print "Will perform stats on all molecule files in the folder  " + fn
print "================================================================================"

#===========================================================================================================================================================================
#Script
#===========================================================================================================================================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Folder structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Try to make result folder
try:
	os.system("mkdir -v "+ fn+fnOut)
except OSError:
	pass
	
#Try to make figs folder
try:
	os.system("mkdir -v "+ fn+fnOut+'figs/')
except OSError:
	pass

#Try to make temp folder
try:
	os.system("mkdir -v "+ fn+fnOut+'figs/temp/')
except OSError:
	pass

#Empty temp folder
os.system('rm -vr '+ fn+fnOut+'figs/temp/*')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
#Get all molecule files
molFiles=pyfrp_misc_module.getSortedFileList(fn,".mol")

#Loop through all molecule files
for f in molFiles:
	
	#Load molecule
	mol=pyfrp_IO_module.loadFromPickle(fn+f)
	
	#Try to make mol temp folder
	try:
		os.system("mkdir -v "+ fn+fnOut+'figs/temp/'+mol.name+"/")
	except OSError:
		pass
	
	#Loop through embryos
	for i,emb in enumerate(mol.embryos):
		
		print "Processing embryo ", emb.name
	
		#Make figure for results
		fig,axes=pyfrp_plot_module.makeSubplot([int(np.ceil(len(emb.fits)/2.)),2],tight=True)
		
		#Plot fits
		for j,fit in enumerate(emb.fits):
			print len(axes), i*2+j
			fit.plotFit(ax=axes[i*2+j])
		
		#Save figure
		fig.savefig(fn+fnOut+'figs/temp/'+mol.name+'/'+emb.name+'.pdf')
		
		
	
		
		plt.close('all')
	
	
	os.system("pdftk "+fn_out+"figs/temp/*pdf cat output "+ fn_out+"figs/"+name+"_allfits.pdf")
	os.system("rm -v "+fn_out+"figs/temp/*pdf")

		
	print "Done with molecule "+ f 

	
	
	
	
print "Done"	
	

