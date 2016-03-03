### Script to sort Theresa files

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import os
import time
from optparse import OptionParser

#Numpy/Scipy
from numpy import *



#===========================================================================================================================================================================
#Defining parsing possibilities
#===========================================================================================================================================================================

usage = "usage: %prog [options]"
parser = OptionParser(usage)
     
parser.add_option("-n", "--name",action="store", type='string', dest="n", help="Name of experiment")
parser.add_option("-d", "--date",action="store", type='string', dest="d", help="Date of experiment")
parser.add_option("-f", "--folder",action="store", type='string', dest="f", help="Data folder")
parser.add_option("-p", "--print",action="store", type='int', dest="p", help="Only print, not move?")


parser.set_defaults(d="")
parser.set_defaults(f="")
parser.set_defaults(n="")
parser.set_defaults(p=0)


#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.f)
date=str(options.d)
name=str(options.n)
onlyprint=bool(int(options.p))

print "================================================================================"
print "Will sort the following folder" + fn
print "name: ", name
print "date: ", date 
print "only print: ", onlyprint
print "================================================================================"


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
	os.mkdir(fn_dest+"/nobeads")
except OSError:	
	print "Could not make folders:"
	print fn_dest+"/nobeads"
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Move and rename embryo folders name
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 2: Move and rename embryo folders name
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

#Get all files in molecule folder
files=os.listdir(fn)

#Get number of experiments in folder
NExp=0
for f in files:
	if "exp" in f and os.path.isdir(fn+"/"+f):
		NExp=NExp+1

#Loop through files
for f in files:
	print
	#Check if folder
	if os.path.isdir(fn+"/"+f):
		
		#Check if exp folder
		if "exp" in f:
			
			#Get number of exp
			nExp=f.strip("exp")
			
			#Fill up number string
			nExp=(len(str(NExp))-len(nExp))*'0'+nExp
			
			#Move lsm files and rename them
			lsmFolder=fn+"/"+f+"/lsm"
			
			try:
				os.mkdir(fn_dest+"/lsm")
			except OSError:	
				print "Could not make folders:"
				print fn_dest+"/lsm"
			
			newlsmFolder=fn_dest+"/lsm/lsm_"+date+"_"+nExp
			
			if onlyprint:
				print "mv -v " + lsmFolder+" "+newlsmFolder
			else:
				os.system("mv -v " + lsmFolder+" "+newlsmFolder)
			
			#New folder name
			new_fn=fn_dest+"/nobeads/"+date+"_"+nExp
			
			#Move files
			if onlyprint:
				print "mv -v " + fn+"/"+f+" "+new_fn
			else:
				os.system("mv -v " + fn+"/"+f+" "+new_fn)
			
			#Rename post folder to recover
			if onlyprint:
				print "mv -v " + new_fn+"/post"+" "+new_fn+"/recover"
			else:
				os.system("mv -v " + new_fn+"/post"+" "+new_fn+"/recover")
			
			#Delete result folder
			if onlyprint:
				print "rm -v -r " + new_fn+"/results"
			else:
				os.system("rm -v -r " + new_fn+"/results")

			
#Move info.txt if existent
if onlyprint:
	print "mv -v " + fn+"/info.txt"+" " +fn_dest+"/info_"+date+".txt"
else:
	os.system("mv -v " + fn+"/info.txt"+" " +fn_dest+"/info_"+date+".txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 5: Clean up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 5: Clean Up
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

if onlyprint:
	print "rm -v -r " + fn_dest+"/results"
else:
	os.system("rm -v -r " + fn_dest+"/results")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 6: Rename main folder
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

print """
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Step 6: Rename main folder
~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

if onlyprint:
	print "mv -v "+fn_dest+" " + root+"/"+name
else:
	os.system("mv -v "+fn_dest+" " + root+"/"+name)


			
			
			
	
	
	