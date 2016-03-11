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
import matplotlib

#===========================================================================================================================================================================
#Defining parsing possibilities
#===========================================================================================================================================================================

usage = "usage: %prog [options]"
parser = OptionParser(usage)
     
parser.add_option("-l", "--file",action="store", type='string', dest="fn", help="Input file")
parser.add_option("-O", "--fnout",action="store", type='string', dest="O", help="Folder for output")
parser.add_option("-N", "--name",type='str',  dest="N", help="Name?")

parser.set_defaults(fn=None)

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
fn_out=str(options.O)
name=str(options.N)


print "================================================================================"
print "Will summarize all datasets specified in file " + fn
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

#---------------------------
#Bookkeeping lists
Ds=[]
Dsterr=[]
names=[]

#---------------------------
#Loop through all molecules specified
for k,fn_mol in enumerate(mols):
	
	print "======================================================================================================================================="
	print "Processing molecule " + fn_mol
	print "======================================================================================================================================="
	
	#Creating molecule object and load molecule file
	mol=molecule("blub")
	mol=mol.load_molecule(fn_mol)
	
	Ds.append(mol.D_mu_av)
	print Ds
	Dsterr.append(mol.D_mu_sterr)
	names.append(mol.name.replace("_100pg","").replace("_"," \n "))
	print names

names=["GFP","secGFP"]

#---------------------------
#Plot Settings
from mpl_toolkits.axes_grid.axislines import SubplotZero

fig_width_pt = 307.28
inches_per_pt = 1.0/72.27
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
	'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
         'xtick.labelsize': 10,
          'ytick.labelsize': 10,
         'text.usetex': True,
         'font.family': 'sans-serif',
         'font.sans-serif': 'Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif',
         'figure.figsize': fig_size}
plt.rcParams.update(params)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Bitstream Vera Sans', 'Lucida Grande', 'Verdana', 'Geneva']  
plt.rcParams['ytick.direction'] = 'out' 

plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!  
]

#---------------------------
#Create figure

fig=plt.figure()
ax=fig.add_subplot(111)

fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(left=0.2)
fig.subplots_adjust(right=0.8)

	
N=shape(Ds)[0]
ind=2.*arange(N)
print ind
width=0.2*2.

Ds=asarray(Ds)
print Ds
Dsterr=asarray(Dsterr)

bar_D=ax.bar(ind-width, Ds, width, color='0.75',yerr=Dsterr,error_kw=dict(ecolor='k'))

ax.set_ylabel("diffusivity (\SI{}{\micro\meter^2/s})")

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_xticks(ind-0.5*width)
ax.set_xticklabels(names)

ax.set_xlim([min(ind)-0.5,max(ind)+0.5])

plt.draw()
plt.savefig('/home/alex_loc/Documents/Talks/PhDSymposium_2015/figs/bar_'+name+'.eps')
plt.show()

plt.draw()
raw_input()
