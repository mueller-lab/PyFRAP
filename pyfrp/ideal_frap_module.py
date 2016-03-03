#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#PyFRAP
import pyfrp_molecule
import pyfrp_IO_module
import pyfrp_plot_module

#Numpy
import numpy as np

#Misc
import time
import os
import copy as cpy

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#resets Analysis options in all different combinations and saves a molecule file for each combination

def resetAnalysisOptionsOfMol(mol,normOpt,flattenOpt,bkgdOpt,gaussianOpt,medianOpt,fn):
	
	#Loop through all options
	for norm in normOpt:
		for bkgd in bkgdOpt:
			for gaussian in gaussianOpt:
				for median in medianOpt:
					
					#Copy embryo
					newMol=cpy.deepcopy(mol)
					
					#Loop through all embryos
					for emb in newMol.embryos:
						
						#Reset analysis options
						emb.analysis.setNorm(bool(norm))
						emb.analysis.setBkgd(bool(bkgd))
						emb.analysis.setFlatten(False)
						emb.analysis.setMedian(bool(median))
						emb.analysis.setGaussian(bool(gaussian))
					
					#Define molecule name
					name,dump=mol.name.split("_n0")
					newMol.name=name+"_n"+str(int(emb.analysis.normOn()))+"_f"+str(int(emb.analysis.flattenOn()))+"_b"+str(int(emb.analysis.bkgdOn()))+"_g"+str(int(emb.analysis.gaussianOn()))+"_m"+str(int(emb.analysis.medianOn()))
	
					#Save
					newMol.save(fn+newMol.name+'.mol')
	
	#Do the same with flattenOptions
	for flatten in flattenOpt:
		for bkgd in bkgdOpt:
			for gaussian in gaussianOpt:
				for median in medianOpt:
					
					#Copy embryo
					newMol=cpy.deepcopy(mol)
					
					#Loop through all embryos
					for emb in newMol.embryos:
						
						#Reset analysis options
						emb.analysis.setNorm(False)
						emb.analysis.setBkgd(bool(bkgd))
						emb.analysis.setFlatten(bool(flatten))
						emb.analysis.setMedian(bool(median))
						emb.analysis.setGaussian(bool(gaussian))
					
					#Define molecule name
					name,dump=mol.name.split("_n0")
					newMol.name=name+"_n"+str(int(emb.analysis.normOn()))+"_f"+str(int(emb.analysis.flattenOn()))+"_b"+str(int(emb.analysis.bkgdOn()))+"_g"+str(int(emb.analysis.gaussianOn()))+"_m"+str(int(emb.analysis.medianOn()))
	
					#Save
					newMol.save(fn+newMol.name+'.mol')	

	return True

#def genEmbryo(name,fnData,fnPre,imagingRadius,sideLengthBleachedMu)

	##Create new embryo
	#emb=mol.newEmbryo(name)
	
	##Set filepath to test dataset
	#emb.setDataFolder(fn+expName+"/"+exp+"/recover/")
	#emb.setPreimage(fn+expName+"/"+exp+"/pre/")
	
	##Update File list
	#emb.updateFileList()
	
	##Set data resolution
	#emb.setDataResMu(dataResMu)

	##Get center and radius of experiment
	#boundarySelector=pyfrp_plot_module.FRAPBoundarySelector(emb)
	#center,imagingRadius=boundarySelector.getResults()
	
	##Define Sidelength of bleached area
	#emb.setSideLengthBleachedMu(sideLengthBleachedMu)
	
	##Define depth of imaging
	#emb.setSliceDepthMu(imagingDepth)
	
	##Define geometry
	#if exp=="embryo":
		#emb.setGeometry2ZebraFishDomeStage(center,imagingRadius)
	#else:
		#emb.setGeometry2Cone(center,imagingRadius,cylinderHeight)
	
	##Update everything in geofile
	#emb.geometry.updateGeoFile()

	##Create default ROIs
	#emb.genDefaultROIs(emb.geometry.getCenter(),emb.geometry.getRadius())

	##Create optimal All ROI
	#emb.getOptimalAllROI()

	##Set frameInterval
	#emb.setFrameInterval(frameInterval)
	
	##Create Simulation
	#sim=emb.newSimulation()

	##Set mesh volSize
	#sim.mesh.setVolSizePx(volSizePx)

	##Generate Mesh
	#sim.mesh.genMesh()
	
	##Compute ROI Idxs
	##emb.computeROIIdxs() #Don't do this here, takes up too much time. Can do this later on a cluster
	
	##Create analysis
	#emb.newAnalysis()

	##Set analysis options
	#emb.analysis.setNorm(norm)
	#emb.analysis.setGaussian(False)
	#emb.analysis.setMedian(True)
	#emb.analysis.setQuad(False)