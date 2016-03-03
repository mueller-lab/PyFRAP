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
import time


from mpl_toolkits.axes_grid.axislines import SubplotZero

fig_width_pt = 307.28
inches_per_pt = 1.0/72.27
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
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
  

plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\usepackage[latin1]{inputenc}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!  
]


fig=plt.figure()
ax=fig.add_subplot(111)

fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(right=0.8)
	
Ds=[11.1273566667,11.212]
prods=[1.77,1.73]
err_Ds=[0.5511440136,0.541]
err_prods=[0.10913,0.1123]

names=["Mueller (2012)*", "PyFRAP"]

N=2
#ind= arange(N)
ind=asarray([0,2])
width=0.3

bar_D=ax.bar(ind-width, Ds, width, color='0.75',yerr=err_Ds,error_kw=dict(ecolor='k'))
ax2=ax.twinx()

bar_prod=ax2.bar(ind, prods, width, color='0.5',yerr=err_prods,error_kw=dict(ecolor='k'))

lg=plt.legend([bar_D,bar_prod],["diffusivity","production"],loc=9,bbox_to_anchor=(0.5, 1))
lg.draw_frame(False)

ax.set_xticks(ind)
ax.set_xticklabels(names)

#ax.set_ylabel("diffusivity ($\mu\mathrm{m}^2/\mathrm{s}$)")
ax.set_ylabel("diffusivity (\SI{}{\micro\meter^2/s})")
ax2.set_ylabel("production ($10^{-4}\, 1/\mathrm{s}$)")

plt.draw()
plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/bar_mueller.eps')
plt.show()