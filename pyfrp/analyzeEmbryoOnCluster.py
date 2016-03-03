#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

import pyfrp_IO_module
import sys
from pyfrp_term_module import *

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

#Pass filename
fn=sys.argv[1]

#Load Embryo
emb=pyfrp_IO_module.loadFromPickle(fn)

#Make sure to have a dense enough mesh
squROI=emb.getROIByName('Bleached Square')
emb.simulation.mesh.setVolSizePx(45.)
emb.simulation.mesh.forceMinMeshDensityInROI(squROI,0.001,stepPercentage=0.1,debug=False,findIdxs=True,method='volSize')

#Image analysis
emb.analysis.run(showProgress=False)

#Simulation
emb.simulation.getOptTvecSim(150.)
emb.simulation.run(showProgress=False)

#Pin Concentrations
bkgdVal,normVal,bkgdValSim,normValSim=emb.computeIdealFRAPPinVals(debug=True,useMin=False,useMax=False,switchThresh=0.95)
emb.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=False)

#Run all fits
for fit in emb.fits:
	fit.run(debug=False)
	
#Save Embryo
emb.save(fn)

#raw_input()	