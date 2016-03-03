###Final script to create molecule files for IdealFRAP Project

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#PyFRAP
import pyfrp_molecule
import pyfrp_IO_module
import pyfrp_plot_module
import ideal_frap_module

#Numpy
import numpy as np

#Misc
import time
import os
from optparse import OptionParser

#===========================================================================================================================================================================
#Defining parsing possibilities
#===========================================================================================================================================================================

usage = "usage: %prog [options]"
parser = OptionParser(usage)
     
parser.add_option("-f", "--fn",action="store", type='string', dest="fn", help="Folder containing molecule data (eg.FITC-dextran-10kDa)")
parser.add_option("-n", "--name",action="store", type='string', dest="name", help="Name of molecule. Equals folder name if not specified.")
parser.add_option("-d", "--depth",type='float', dest="depth", help="Imaging depth. Default=40.")
parser.add_option("-r", "--res",type='float', dest="res", help="Resolution in mu. Default=566.79")
parser.add_option("-i", "--dt", type='float', dest="dt", help="Frame interval in s. Default=1.")
parser.add_option("-v", "--vol",type='float',  dest="vol", help="Mesh cell volume?  Default=23")
parser.add_option("-t", "--steps",type='int',  dest="steps", help="Simulation steps?  Default=3000")
parser.add_option("-s", "--sl",type='float',  dest="sl", help="Sidelength of bleached region in mu?  Default=192.68")
parser.add_option("-c", "--copy",type='int',  dest="copy", help="Only copy?  Default=False")
parser.add_option("-C", "--nocopy",type='int',  dest="nocopy", help="Only generate?  Default=False")

#Set defaults 
parser.set_defaults(fn=None)
parser.set_defaults(name=None)
parser.set_defaults(depth=40.)
parser.set_defaults(res=566.79)
parser.set_defaults(dt=1.)
parser.set_defaults(vol=23)
parser.set_defaults(norm=0)
parser.set_defaults(steps=3000)
parser.set_defaults(sl=192.68)
parser.set_defaults(copy=0)
parser.set_defaults(nocopy=0)

#Getting parameters
(options, args) = parser.parse_args()     

fn=str(options.fn)
norm=bool(options.norm)
onlyCopy=bool(options.copy)
stepsSim=int(options.steps)
volSizePx=float(options.vol)
dataResMu=float(options.res)
frameInterval=float(options.dt)
imagingDepth=float(options.depth)
name=str(options.name)
sideLengthBleachedMu=float(options.sl)
noCopy=bool(options.nocopy)

print "================================================================================"
print "Will create molecule files in the folder  " + fn
print "Will add the following parameters:"

print "name: ", name
print "onlyCopy: ", onlyCopy
print "noCopy: ", noCopy

print "Data:"
print "imagingDepth: ", imagingDepth
print "dataResMu: ", dataResMu
print "sideLengthBleachedMu: ", sideLengthBleachedMu
print "frameInterval: ", frameInterval

print "Analysis:"
print "norm: ", norm

print "Simulation"
print "volSizePx: ", volSizePx
print "stepsSim: ", stepsSim

print "================================================================================"

#===========================================================================================================================================================================
#Default parameters
#===========================================================================================================================================================================

cylinderHeight=90.3332804037
coneUpperRadius=635.3/2.
coneLowerRadius=448.5/2.
molecules=[]

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

#Start timer
startInit=time.clock()

#Check for molecule name
if name=='None':
	if fn[-1]=="/":
		root,fnName=os.path.split(fn[:-1])
	else:
		root,fnName=os.path.split(fn)
		fn=fn+"/"
		
	name=fnName
	if fnName=='':
		print "fnName empty"
		raw_input()

#Define different experiment types
expNames=["nobeads","beads","embryo"]

#Loop through all experiment types, if experiment type exists create molecule
for typeCount,expName in enumerate(expNames):
	
	print
	
	#Check if this type of experiment was done
	if not os.path.isdir(fn+expName):
		continue
	
	#Check if we need to create master molecule file
	if not onlyCopy:
		
		#Create molecule object
		mol=pyfrp_molecule.molecule(name+"_"+expName)
		
		#Check for all existent experiments
		experiments=os.listdir(fn+expName)
		
		#Create embryo object for each embryo and set options
		for expCount,exp in enumerate(experiments):
			
			if expCount>2:
				break
			
			#Create new embryo
			emb=mol.newEmbryo(exp)
			
			#Set filepath to test dataset
			emb.setDataFolder(fn+expName+"/"+exp+"/recover/")
			
			#Update File list
			emb.updateFileList()
			
			#Set data resolution
			emb.setDataResMu(dataResMu)

			#Get center and radius of experiment
			boundarySelector=pyfrp_plot_module.FRAPBoundarySelector(emb)
			center,imagingRadius=boundarySelector.getResults()
			
			#Define Sidelength of bleached area
			emb.setSideLengthBleachedMu(sideLengthBleachedMu)
			
			#Define depth of imaging
			emb.setSliceDepthMu(imagingDepth)
			
			#Define geometry
			if exp=="embryo":
				emb.setGeometry2ZebraFishDomeStage(center,imagingRadius)
			else:
				emb.setGeometry2Cone(center,coneUpperRadius,coneLowerRadius,-cylinderHeight)
				
				#In case of cone, return the slice height from cone geometry
				emb.sliceHeightPx=emb.geometry.computeSliceHeightFromRadius(imagingRadius)
				emb.sliceDepthPx=-emb.sliceHeightPx
				emb.sliceDepthMu=emb.sliceDepthPx*emb.convFact
			
			#Update everything in geofile
			emb.geometry.updateGeoFile()

			#Create default ROIs
			emb.genDefaultROIs(emb.geometry.getCenter(),imagingRadius)
		
			#Create optimal All ROI
			emb.getOptimalAllROI()

			#Set frameInterval
			emb.setFrameInterval(frameInterval)
			
			#Create Simulation
			sim=emb.newSimulation()

			#Set mesh volSize
			sim.mesh.setVolSizePx(volSizePx)
			
			#Set slice ROI to be refined
			sl=emb.getROIByName("Slice")
			sl.refineInMesh(debug=False,factor=5.,addZ=5.,findIdxs=False)
			
			#Generate Mesh
			sim.mesh.genMesh()
			
			#Compute ROI Idxs
			#emb.computeROIIdxs() #Don't do this here, takes up too much time. Can do this later on a cluster
			
			#Create analysis
			emb.newAnalysis()

			#Set analysis options
			emb.analysis.setNorm(False)
			emb.analysis.setGaussian(False)
			emb.analysis.setMedian(False)
			emb.analysis.setQuad(False)
			emb.analysis.setBkgd(False)
			emb.analysis.setFlatten(False)
			emb.analysis.setFnBkgd("../Ideal_FRAP/backgroundData/20160105_bkgd/01/")
			emb.analysis.setFnFlatten("../Ideal_FRAP/flatteningData/20160105_flat/01/")
			emb.analysis.setFnPre(fn+expName+"/"+exp+"/pre/")
			
			#Append fits
			for fitDegr in [0,1]:
				for fitProd in [0,1]:
					for eq in [0,1]: 
					
						fit=emb.newFit('fit'+'_p'+str(fitProd)+'_d'+str(fitDegr)+'_eq'+str(eq))
						fit.addROIByName('Bleached Square')
						fit.addROIByName('Slice')
						
						fit.fitDegr=bool(fitDegr)
						fit.fitProd=bool(fitProd)
						
						fit.equOn=bool(eq)
						fit.fitPinned=bool(eq)
						
						fit.setOptMeth('Constrained Nelder-Mead')
				
			print "Created embryo " , emb.name
		
		#Include analysis options in molecule name
		mol.name=mol.name+"_n"+str(int(emb.analysis.normOn()))+"_f"+str(int(emb.analysis.flattenOn()))+"_b"+str(int(emb.analysis.bkgdOn()))+"_g"+str(int(emb.analysis.gaussianOn()))+"_m"+str(int(emb.analysis.medianOn()))
		
		print "Initialized molecule ", mol.name, "  after ", time.clock()-startInit
		
		#Save
		mol.save(fn+mol.name+'.mol')
	
	else:
		#Load
		mol=pyfrp_IO_module.loadFromPickle(fn+name+"_"+expName+"_n"+str(0)+"_f"+str(0)+"_b"+str(0)+"_g"+str(0)+"_m"+str(0)+".mol")
	
	#Make all other molecule files
	if not noCopy:
		ideal_frap_module.resetAnalysisOptionsOfMol(mol,[0,1],[0,1],[0,1],[0,1],[0,1],fn)
	
print "Done"	
	

