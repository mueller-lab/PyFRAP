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
     
parser.add_option("-l", "--file",action="store", type='string', dest="fn", help="Input file")
parser.add_option("-n", "--nfit",action="store", type='string', dest="n", help="Selected fit")
parser.add_option("-O", "--fnout",action="store", type='string', dest="O", help="Folder for output")
parser.add_option("-e", "--embs",type='str',  dest="e", help="embryos?")
parser.add_option("-N", "--name",type='str',  dest="N", help="Name?")
parser.add_option("-x", "--exclude",type='str',  dest="x", help="exclude?")


parser.set_defaults(fn=None)

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
embs=str(options.e)
nfit=int(options.n)
fn_out=str(options.O)
name=str(options.N)
exclude=str(options.x)


print "================================================================================"
print "Will summarize all datasets specified in file " + fn
print "Will do the following:"
print "Embryos: ", embs
print "Exclude: ", exclude
print "nfit: ", nfit
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
#Parameters
#===========================================================================================================================================================================

#Read in embryos to analyze
if embs=='None':
	embs=len(mols)*"[],"
	embs="["+embs[:-1]+"]"
	
embs,dump=pyfrp_misc.str2list(embs)

#Read in embryos to analyze
if exclude=='None':
	exclude=len(mols)*"[],"
	exclude="["+exclude[:-1]+"]"
	
exclude,dump=pyfrp_misc.str2list(exclude)

#===========================================================================================================================================================================
#Bookkeeping lists
#===========================================================================================================================================================================

names=[]
D=[]
Rsq=[]
ssd=[]
Dstd=[]
Dsterr=[]
allrows=[]

allrows_emb=[]

allrows_emb.append(["name","D_opt_mu","Rsq","ssd","fit_pinned","equ_on","fit_squ","fit_out","fit_slice"])

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

#Loop through all molecules specified
for k,fn_mol in enumerate(mols):
	
	print "======================================================================================================================================="
	print "Processing molecule " + fn_mol
	print "======================================================================================================================================="
	
	#Creating molecule object and load molecule file
	mol=molecule("blub")
	mol=mol.load_molecule(fn_mol)
	
	#Take all embryos if embs is empty
	if len(embs[k])==0:
		embs[k]=range(len(mol.embryos))
	
	#Empty selected fits
	mol.sel_fits=[]
	
	#Append name to allrows
	allrows_emb.append([mol.name])
	
	#Loop through embryos
	for n,emb in enumerate(mol.embryos):
		
		#Check if we are supposed to analyze this embryo
		if n in embs[k] and n not in exclude[k]:
		
			#Check if it is analyzed
			if not len(emb.fits)>0:
				print "Warning, embryo " + emb.name + " of molecule " + fn_mol + "  is not analyzed yet"
				break
			
			#Grab fit
			fit=emb.fits[nfit]
			
			#Add fit to sel_fits
			mol.sel_fits.append(fit)
			
			#Save figures
			if fn_out[-1]=="/":
				pass
			else:
				fn_out=fn_out+"/"
			
			fn_save=fn_out+"figs/temp/"+mol.name+"_"+emb.name+"_norm"+str(int(emb.norm_by_pre))+"_fit"+str(nfit)+".pdf"
				
			if fit.fit_pinned:
				fit.save_plot_fit_pinned(fn_save,title=mol.name+"_"+emb.name+"_norm"+str(int(emb.norm_by_pre))+"_fit"+str(nfit))
			else:	
				fit.save_plot_fit_unpinned(fn_save,title=mol.name+"_"+emb.name+"_norm"+str(int(emb.norm_by_pre))+"_fit"+str(nfit))
			
			
			
			allrows_emb.append([emb.name,fit.D_opt_mu,fit.Rsq,fit.ssd,fit.fit_pinned,fit.equ_on,fit.fit_squ,fit.fit_out,fit.fit_slice])
			
	mol.sumup_results()
	
	names.append(mol.name)
	D.append(mol.D_mu_av)
	Rsq.append(mol.Rsq_av)
	ssd.append(mol.ssd_av)
	Dstd.append(mol.D_mu_std)
	Dsterr.append(mol.D_mu_sterr)
	print mol.prod_av
	raw_input()
	wrow=[mol.name,mol.D_mu_av,mol.D_mu_std,mol.D_mu_sterr,mol.Rsq_av,mol.ssd_av]
	
	allrows.append(wrow)
	
	mol.save_molecule(fn_mol)
	
#Concatenate pdfs
os.system("pdftk "+fn_out+"figs/temp/*pdf cat output "+ fn_out+"figs/"+name+"_allfits.pdf")
os.system("rm -v "+fn_out+"figs/temp/*pdf")

#Write molecule averages to csv sheet	
fn_csv=fn_out+"/"+name+"_summary_molecules.csv"
wfile=csv.writer(open(fn_csv,'wb'), delimiter=';')

nrow=["name","D_mu_av","D_mu_std","D_mu_sterr","Rsq_av","ssd_av"]
	
wfile.writerow(nrow)

for r in allrows:
	wfile.writerow(r)

#Write single embryo results to csv sheet	
fn_csv=fn_out+"/"+name+"_summary_embryos.csv"
wfile=csv.writer(open(fn_csv,'wb'), delimiter=';')

for r in allrows_emb:
	wfile.writerow(r)

print 	
print "Done"
	
