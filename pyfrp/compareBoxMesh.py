#Final test script for PyFRAP
import pyfrp_molecule
import pyfrp_IO_module
import pyfrp_plot_module

import numpy as np

import time

load=1
runNoBox=0
runBox=1

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
	#emb.save(fn="../TestDataset/embryoFiles/TestEmbryo.emb")
	
#NoBox
if runNoBox:	
	
	#Analysis
	emb.analysis.run(debug=False)	
	mol.save('../TestDataset/'+mol.name+'_nobox.mol')
	mol.save('../TestDataset/'+mol.name+'.mol')
	
	#Simulation
	sim=emb.simulation
	sim.setTimesteps(3000)
	sim.setD(10.)
	sim.getOptTvecSim(150.)
	
	sim.run()
	mol.save('../TestDataset/'+mol.name+'_nobox.mol')
	
	#Pin
	bkgdVal,normVal,bkgdValSim,normValSim=emb.computeIdealFRAPPinVals(debug=True,useMin=False,useMax=False,switchThresh=0.95)
	emb.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=False)
	
	#Run fit
	emb.fits[0].run(debug=False)
	emb.fits[0].printResults()
	mol.save('../TestDataset/'+mol.name+'_nobox.mol')

#Box
if runBox:
	
	#Refine and run again
	sl=emb.getROIByName("Slice")
	al=emb.getROIByName("All")
	sl.refineInMesh(debug=True,factor=5.,addZ=5.)
	emb.computeROIIdxs()
	#print len(sl.meshIdx)/float(len(al.meshIdx)), abs(sl.getZExtend()[0]-sl.getZExtend()[1])/90.3

	#Simulation
	sim=emb.simulation
	sim.setTimesteps(3000)
	sim.getOptTvecSim(150.)
	sim.toLogTimeScale()
	sim.setD(10.)

	sim.run()
	mol.save('../TestDataset/'+mol.name+'_box.mol')
	
	#Pin
	bkgdVal,normVal,bkgdValSim,normValSim=emb.computeIdealFRAPPinVals(debug=True,useMin=False,useMax=False,switchThresh=0.95)
	emb.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=False)
	
	#Run fit
	emb.fits[0].run(debug=False)
	emb.fits[0].printResults()
	mol.save('../TestDataset/'+mol.name+'_box.mol')
	
molBox=pyfrp_IO_module.loadFromPickle('../TestDataset/TestMolecule_box.mol')
molNoBox=pyfrp_IO_module.loadFromPickle('../TestDataset/TestMolecule_nobox.mol')

embBox=molBox.embryos[0]
embNoBox=molNoBox.embryos[0]

squBox=embBox.getROIByName('Bleached Square')
squNoBox=embNoBox.getROIByName('Bleached Square')

sliceBox=embBox.getROIByName('Slice')
sliceNoBox=embNoBox.getROIByName('Slice')

allBox=embBox.getROIByName('All')
allNoBox=embNoBox.getROIByName('All')

print "Interpolation Error squBox:", squBox.getInterpolationError()
print "Interpolation Error squNoBox:", squNoBox.getInterpolationError()

print "Interpolation Error sliceBox:", sliceBox.getInterpolationError()
print "Interpolation Error sliceNoBox:", sliceNoBox.getInterpolationError()

fitBox=embBox.fits[0]
fitNoBox=embNoBox.fits[0]

#sim=embBox.simulation

#sim.showInterpolatedIC()


#ax=embNoBox.plotAllData()
#ax=embNoBox.plotAllSim(ax=ax)

#ax=embBox.plotAllDataPinned()
#ax=embBox.plotAllSimPinned(ax=ax)

#embBox.simulation.compareICInterpolation(roi=allBox)

#allBox.plotSolutionVariable(sim.IC,nlevels=100)

#ax=embBox.plotAllData()
#ax=embBox.plotAllSim(ax=ax)

#ax=embBox.plotAllDataPinned()
#ax=embBox.plotAllSimPinned(ax=ax)

fitBox.plotFit()
fitNoBox.plotFit()



raw_input()


	