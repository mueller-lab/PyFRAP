from fipy import *
import math 
import os
from numpy import *
from pyfrp_sim_module import *
from pyfrp_misc_module import *
from pyfrp_img_module import *
from pyfrp_fit_module import *
from pyfrp_stats_module import *
from pyfrp_gmsh_module import *

from embryo import *
from molecule import *

import time

datasets=["A488-dextran-3kDa_20140612","A488-dextran-3kDa_1uM_agarose_1Percent_20140616","A488-dextran-10kDa_20140610","A488-dextran-10kDa_agarose_1percent_20140616"]

Ds_FCS=[141,107,80,64]
Ds_FCS=asarray(Ds_FCS)

Ds=[]

err_Ds=[]
for dataset in datasets:
	#Loading molecule
	mol=molecule(dataset)
	mol=mol.load_molecule("/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/"+dataset+"/results/"+dataset+".pk")
	
	
	Ds_temp=[]
	
	
	mol.sel_fits=[]
			
	for emb in mol.embryos:
		mol.sel_fits.append(emb.fits[0])	
		Ds_temp.append(emb.fits[0].D_opt_mu)
	
	#err_D.append(std(Ds_temp)/sqrt(len(Ds_temp)))
	
	err_Ds.append(std(Ds_temp))
	mol.sumup_results()
		
	Ds.append(mol.D_mu_av)	
	


from mpl_toolkits.axes_grid.axislines import SubplotZero

fig_width_pt = 307.28
inches_per_pt = 1.0/72.27
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_height = 2.02
fig_size =  [fig_width,fig_height]
print fig_size
params = {'backend': 'ps',
	'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
         'xtick.labelsize': 10,
          'ytick.labelsize': 10,
         'text.usetex': True,
         'font.family': 'sans-serif',
         'font.sans-serif': 'Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif',
         'figure.figsize': fig_size}
plt.rcParams.update(params)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Bitstream Vera Sans', 'Lucida Grande', 'Verdana', 'Geneva']  
plt.rcParams['ytick.direction'] = 'out' 

plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!  
]


fig=plt.figure()
ax=fig.add_subplot(111)

fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(left=0.2)
fig.subplots_adjust(right=0.8)


print Ds[0]/Ds[2], Ds[1]/Ds[3]

Ds=asarray(Ds)

err_Ds=asarray(err_Ds)


names=["3 kDa","3 kDa \n in agarose","10 kDa","10 kDa \n in agarose"]

N=4
ind= arange(N)
#ind=asarray([0,2])
width=0.2

bar_D=ax.bar(ind-width, Ds, width, color='0.75',yerr=err_Ds,error_kw=dict(ecolor='k'))

#bar_FCS=ax.bar(ind, Ds_FCS, width, color='c')

#lg=plt.legend([bar_D,bar_prod],["diffusivity","production"],loc=9,bbox_to_anchor=(0.5, 1))
#lg.draw_frame(False)

ax.set_xticks(ind-0.5*width)
ax.set_xticklabels(names)

ax.set_ylabel("diffusivity (\SI{}{\micro\meter^2/s})")

ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.draw()
plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/bar_vitro.eps')
plt.show()