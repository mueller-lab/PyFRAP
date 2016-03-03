#Final test script for PyFRAP
import pyfrp_molecule
import pyfrp_IO_module
import pyfrp_plot_module

import numpy as np

import time

load=1
analyze=0
simulate=0
fitemb=1
loadOld=1
loadOldIntoEmb=1
pinemb=0

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
	emb.setDataFolder("../Ideal_FRAP/Gary/Fscin40kDa_500nM/nobeads/20150818_01/recover/")
	emb.setPreimage("../Ideal_FRAP/Gary/Fscin40kDa_500nM/nobeads/20150818_01/pre/Fscin40kDa_500nM_01_pre_t000.tif")
	
	#Update File list
	emb.updateFileList()
	
	#Set data resolution
	emb.setDataResMu(dataResMu)
	
	#Define Sidelength of bleached area
	emb.setSideLengthBleachedMu(192.68488493739596)
	
	#Define depth of imaging
	emb.setSliceDepthMu(imagingDepth)
	
	#Define Geometry
	#emb.setGeometry2Cylinder(center,imagingRadius,cylinderHeight)
	emb.setGeometry2CylinderQuad(center,imagingRadius,cylinderHeight)
	
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
	sim.mesh.setVolSizePx(20.)

	#Generate Mesh
	sim.mesh.genMesh()
		
	#Compute ROI Idxs
	emb.computeROIIdxs()
	
	#show Idxs of ROIs	
	#emb.showAllROIIdxs()
	
	emb.makeQuadReducable(auto=True)

	#Create analysis
	emb.newAnalysis()

	#Set analysis options
	emb.analysis.setNorm(False)
	emb.analysis.setGaussian(True)
	emb.analysis.setMedian(False)
	emb.analysis.setQuad(True)
	
	print "Initialized embryo in ", time.clock()-startInit
	
	#Save
	mol.save('testscripts/'+mol.name+'.mol')
	
else:
	
	#Load
	mol=pyfrp_IO_module.loadFromPickle('testscripts/TestMolecule.mol')
	emb=mol.embryos[0]

#Analyze Dataset
if analyze:
	emb.analysis.run(debug=False)	
	mol.save('testscripts/'+mol.name+'.mol')

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
	mol.save('testscripts/'+mol.name+'.mol')

#Load Data for comparison
if loadOld or loadOldIntoEmb:
	tVecDataOld=np.loadtxt('../comparePyFRAPs/tvec_data.txt')
	tVecSimOld=np.loadtxt('../comparePyFRAPs/tvec_sim.txt')

	squDataOld=np.loadtxt('../comparePyFRAPs/squ_data.txt')
	sliceDataOld=np.loadtxt('../comparePyFRAPs/slice_data.txt')
	squSimOld=np.loadtxt('../comparePyFRAPs/squ_sim.txt')
	sliceSimOld=np.loadtxt('../comparePyFRAPs/slice_sim.txt')

	squDataPinnedOld=np.loadtxt('../comparePyFRAPs/squ_data_pinned.txt')
	sliceDataPinnedOld=np.loadtxt('../comparePyFRAPs/slice_data_pinned.txt')
	squSimPinnedOld=np.loadtxt('../comparePyFRAPs/squ_sim_pinned.txt')
	sliceSimPinnedOld=np.loadtxt('../comparePyFRAPs/slice_sim_pinned.txt')

#Load old data into embryo
if loadOldIntoEmb:
	
	emb.tvecData=tVecDataOld
	emb.simulation.tvecSim=tVecSimOld
	
	squROI=emb.getROIByName('Bleached Square')
	sliceROI=emb.getROIByName('Slice')
	
	squROI.dataVec=squDataOld
	sliceROI.dataVec=sliceDataOld
	
	squROI.simVec=squSimOld
	sliceROI.simVec=sliceSimOld
	
	squROI.dataVecPinned=squDataPinnedOld
	sliceROI.dataVecPinned=sliceDataPinnedOld
	
	squROI.simVecPinned=squSimPinnedOld
	sliceROI.simVecPinned=sliceSimPinnedOld

#Pin 
if pinemb:
	#Pin
	bkgdVal,normVal,bkgdValSim,normValSim=emb.computeIdealFRAPPinVals(debug=True,useMin=False,useMax=False,switchThresh=0.95)
	emb.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=False)
	
	ax=squROI.plotDataPinned()
	ax=pyfrp_plot_module.plotTS(tVecDataOld,squDataPinnedOld,color='r',linestyle='--',ax=ax)
	
	ax=squROI.plotSimPinned()
	ax=pyfrp_plot_module.plotTS(tVecSimOld,squSimPinnedOld,color='r',linestyle='--',ax=ax)
	
	raw_input()
	
	
#Fit
if fitemb:
	
	
	#Get ROIs
	squROI=emb.getROIByName('Bleached Square')
	sliceROI=emb.getROIByName('Slice')
	
	#Create fit and add ROIs
	fit=emb.newFit('TestFit')
	fit.addROIByName('Bleached Square')
	fit.addROIByName('Slice')
	
	fit.setOptMeth('Constrained Nelder-Mead')
	
	#Run fit
	fit.run(debug=True)
	
	fit.printResults()
	
	fit.plotFit()

#print squDataOld[0], squROI.dataVec[0]

##Compare Data Series
#ax=squROI.plotData()
#ax=sliceROI.plotData(ax=ax)

#ax=pyfrp_plot_module.plotTS(tVecDataOld,squDataOld,color='r',linestyle='--',ax=ax)
#ax=pyfrp_plot_module.plotTS(tVecDataOld,sliceDataOld,color='r',linestyle='--',ax=ax)

##Compare Simulation Series
#ax=squROI.plotSim(ax=ax)
#ax=sliceROI.plotSim(ax=ax)

#ax=pyfrp_plot_module.plotTS(tVecSimOld,squSimOld,color='r',linestyle='--',ax=ax)
#ax=pyfrp_plot_module.plotTS(tVecSimOld,sliceSimOld,color='r',linestyle='--',ax=ax)






#ax=emb.plotAllSimPinned()
#emb.plotAllDataPinned(ax=ax)

#ax=emb.plotAllSim()
#emb.plotAllData(ax=ax)

#fit.plotFit()
#fit.printResults()

raw_input()



