import pyfrp_IO_module

mol=pyfrp_IO_module.loadFromPickle("../../maria_lft1_FRAP/lft1GFP_morpho_antiGFPmCherry_FRAP_20160203/lft1GFP_morpho_antiGFPmCherry_FRAP_20160203.mol")

fit=mol.embryos[0].fits[1]

fit.x0=[1,0,0]



fit.run(debug=True)

