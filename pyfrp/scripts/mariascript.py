import pyfrp_IO_module

mol=pyfrp_IO_module.loadFromPickle("../../maria_lft1_FRAP/lft1GFP_morpho_antiGFPmCherry_FRAP_20160203/lft1GFP_morpho_antiGFPmCherry_FRAP_20160203.mol")

emb=mol.embryos[0]

emb.simulation.mesh.setVolSizePx(50.)

axes=emb.simulation.mesh.plotDensity()

sl=emb.getROIByName('Slice')

#emb.simulation.mesh.printStats()

sl.refineInMesh(debug=True,factor=5.)

#emb.simulation.mesh.printStats()

emb.simulation.mesh.plotDensity(axes=axes,color='r')

emb.simulation.mesh.plotMesh()

raw_input()

#squ.pinAllTS()
#emb.showDataImg()

