#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Gmsh module for PyFRAP toolbox, including following functions:

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
from numpy import *

#Misc
import os
import shutil
from tempfile import mkstemp
import subprocess
import time
import shlex

#PyFRAP
import pyfrp_gmsh_IO_module
from pyfrp_term_module import *
		
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates cylinder .geo file and remeshes

def updateCylinderGeo(fn,radius,height,center,run=True,debug=False):
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"radius",radius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"height",height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])

	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates cone .geo file and remeshes

def updateConeGeo(fn,upperRadius,lowerRadius,height,center,run=True,debug=False):
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"upper_radius",upperRadius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"lower_radius",lowerRadius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"height",height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])

	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates cylinder .geo file and remeshes

def updateBallGeo(fn,radius,center,run=True,debug=False):
	
	v=5*int(debug)
		
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"radius",radius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])

	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates slice properties in .geo file and remeshes

def updateSliceGeo(fn,sl_height,sl_width,sl_extend,run=True,debug=False):
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"sl_refined",sl_extend)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"sl_height",sl_height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"sl_width",sl_width)
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates refined volume in .geo file and remeshes

def updateRefinedVol(fn,vol,run=True,debug=False):

	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"volSize_fine",vol)
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates cylinder .geo file with refined slice and remeshes

def updateRefinedCylinderGeo(fn,radius,height,center,volSize_px,sl_height,sl_width,sl_extend,vol,run=True,debug=False):
	
	v=5*int(debug)
	
	UpdateCylinderGeo(fn,radius,height,center,volSize_px,run=False,debug=debug)
	UpdateSliceGeo(fn,sl_height,sl_width,sl_extend,run=False,debug=debug)
	updateRefinedVol(fn,vol,run=False,debug=debug)
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates dome .geo file and remeshes

def updateDomeGeo(fn,radius,slice_height,center,run=False,debug=False):
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"radius",radius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"slice_height",slice_height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Updates volSize in  .geo file and remeshes

def updateVolSizeGeo(fn,volSize_px,run=False,debug=False):
	
	v=5*int(debug)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"volSize_px",volSize_px)
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Refines gmsh mesh

def refineMsh(fn,debug=False):
	v=5*int(debug)
	
	gmshBin=getGmshBin()
	os.system(gmshBin+" -v "+ str(v) + " -refine " + fn)
	return fn

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Run gmsh

def runGmsh(fn,fnOut=None,debug=False,redirect=False,fnStout='meshfiles/gmshLogs/gmsh.stout',fnSterr='meshfiles/gmshLogs/gmsh.sterr',volSizeMax=None):
	
	v=5*int(debug)
	
	#Define which command to execute
	gmshBin=getGmshBin()
	cmd = gmshBin + " -v " + str(v) +" -3 -optimize -algo del3d"
	if volSizeMax!=None:
		cmd = cmd + " -clmax " + str(volSizeMax)
	if fnOut!=None:
		cmd = cmd + " -o " + fnOut + " "
	cmd = cmd+ " " + fn
	
	#Print out what will be done
	if debug:
		print "Will execute:"
		print cmd
	
	#Split command in list for subprocess
	args = shlex.split(cmd)
	
	#redirect stdout and stderr if selected
	if redirect:
		stoutFile = open(fnStout,'wb')
		sterrFile = open(fnSterr,'wb')
	else:	
		stoutFile = None
		sterrFile = None
		
	#Call gmsh via subprocess and wait till its done
	try:
		p = subprocess.Popen(args,stdout=stoutFile,stderr=sterrFile)
		p.wait()
	except:
		printError("Gmsh is not running properly, something is wrong.")

	return fn

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gets gmsh executable from path configuration file

def getGmshBin(fnPath="Configurations/paths",identifier="gmshBin"):
	
	if not os.path.isfile(fnPath):
		printWarning(fnPath + " does not exist. Check your paths definition file. Will return 'gmsh'.")
		return "gmsh"
	
	path=None
	
	f = open (fnPath,'rb')
	for line in f:
		if line.strip().startswith(identifier):
			ident,path=line.split('=')
			path=path.strip()
			break
	
	if path==None:
		printWarning("There is no line starting with ", identifier+"= in ", fnPath, ". Will return 'gmsh'.")
		return "gmsh"

	path=os.path.expanduser(path)
	
	return path


	