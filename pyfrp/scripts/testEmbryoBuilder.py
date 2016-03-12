import pyfrp.modules.pyfrp_misc_module as pmm
import os

fn='/mnt/gor_pyfrap/data/Theresa/FITC-dextran-4kDa_sigma_15uM/lsm/'

files=os.listdir(fn)



for f in files:
	if os.path.isdir(fn+f):
		fnDest=fn.replace('lsm/','nobeads/')+f.replace('lsm_','')
		r=pmm.buildEmbryoWizard(fn+f,'lsm','myembryo',fnDest=fnDest,createEmbryo=False)
		if r==-1:
			raw_input()
