#Final test script for PyFRAP
import pyfrp_molecule
import pyfrp_IO_module
import pyfrp_plot_module

import numpy as np

import time

load=0
analyze=1
simulate=1
fitemb=1
pinemb=1


if not load:
	startInit=time.clock()
	
	#Define microscope parameters
	center=[256,256]
	cylinderHeight=90.3332804037
	imagingRadius=294.103867936
	imagingDepth=40
	dataResMu=566.79
	frameInterval=1
	
	#Create test molecule
	mol=pyfrp_molecule.molecule("TestMolecule")
	
	#Create test embryo
	emb=mol.newEmbryo("TestEmbryo")
	
	#Set filepath to test dataset
	emb.setDataFolder("../TestDataset/nobeads/recover/")

	#Update File list
	emb.updateFileList()
	
	#Set data resolution
	emb.setDataResMu(dataResMu)
	
	#Define Sidelength of bleached area
	emb.setSideLengthBleachedMu(192.68488493739596)
	
	#Define depth of imaging
	emb.setSliceDepthMu(imagingDepth)
		
	#Define Geometry
	emb.setGeometry2Cylinder(center,imagingRadius,cylinderHeight)
	
	#Update everything in geofile
	emb.geometry.updateGeoFile()
	
	#Create default ROIs
	emb.genDefaultROIs(emb.geometry.getCenter(),emb.geometry.getRadius())
	
	#Create optimal All ROI
	emb.getOptimalAllROI()
	
	#Set frameInterval
	emb.setFrameInterval(frameInterval)

	#Create Simulation
	sim=emb.newSimulation()

	#Set mesh volSize
	sim.mesh.setVolSizePx(30.)

	#Generate Mesh
	sim.mesh.genMesh()
	
	#Compute ROI Idxs
	emb.computeROIIdxs()
	
	#show Idxs of ROIs	
	#emb.showAllROIIdxs()
	
	#Create analysis
	emb.newAnalysis()
	
	#Set additional data
	emb.analysis.setFnPre("../TestDataset/nobeads/pre/Fscin40kDa_500nM_01_pre_t000.tif")
	emb.analysis.setFnFlatten("../TestDataset/flatten/")
	emb.analysis.setFnBkgd("../TestDataset/bkgd/")
	
	#Set analysis options
	emb.analysis.setNorm(False)
	emb.analysis.setGaussian(False)
	emb.analysis.setMedian(True)
	emb.analysis.setQuad(False)
	emb.analysis.setBkgd(True)
	emb.analysis.setFlatten(True)
	
	emb.analysis.setMedianRadius(1.)
	
	#Get ROIs
	squROI=emb.getROIByName('Bleached Square')
	sliceROI=emb.getROIByName('Slice')
	
	#Create fit and add ROIs
	fit=emb.newFit('TestFit')
	fit.addROIByName('Bleached Square')
	fit.addROIByName('Slice')
	
	fit.setOptMeth('Constrained Nelder-Mead')
	
	print "Initialized embryo in ", time.clock()-startInit
	
	#Save
	mol.save('../TestDataset/'+mol.name+'.mol')
	
else:
	
	#Load
	mol=pyfrp_IO_module.loadFromPickle('../TestDataset/TestMolecule.mol')
	emb=mol.embryos[0]
	emb.save(fn="../TestDataset/embryoFiles/TestEmbryo.emb")

#Analyze Dataset
if analyze:	
	
	emb.analysis.printProcess()
	emb.analysis.run(debug=False)	
	mol.save('../TestDataset/'+mol.name+'.mol')

#Simulate
if simulate:
	
	sim=emb.simulation
	
	#Set Timesteps
	sim.setTimesteps(3000)
	
	#Set diffusion coefficient
	sim.setD(10.)
	
	#Generate Optimal Simulation TimeVector
	sim.getOptTvecSim(150.)
	
	#Set Mesh Properties
	#sim.mesh.setVolSizePx(22.)
	#sim.mesh.genMesh()
	
	sim.run()
	mol.save('../TestDataset/'+mol.name+'.mol')


#Pin 
if pinemb:
	
	#Pin
	bkgdVal,normVal,bkgdValSim,normValSim=emb.computeIdealFRAPPinVals(debug=True,useMin=False,useMax=False,switchThresh=0.95)
	emb.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=False)
	
#Fit
if fitemb:
	
	#Run fit
	fit.run(debug=False)
	
	fit.printResults()
	


