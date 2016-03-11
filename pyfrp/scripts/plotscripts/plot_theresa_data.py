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
import csv

fn_lit="../Ideal_FRAP/Theresa/fitc_dextrans_literature.csv"
fn_alex="../Ideal_FRAP/Theresa/fitc_dextrans_alex.csv"

lit=loadtxt(fn_lit,delimiter=',')
alex=loadtxt(fn_alex,delimiter=',')

#lit=log(lit)
#alex=log(alex)

print shape(lit)

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
plt.rcParams['ytick.direction'] = 'out' 

plt.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!  
]



A = vstack([log(lit[:,0]), ones(len(log(lit[:,0])))]).T
m, c = linalg.lstsq(A, log(lit[:,1]))[0]

fig=plt.figure()
#fig.show()
ax=fig.add_subplot(111)

fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(left=0.2)
fig.subplots_adjust(right=0.8)



#ax.scatter(lit[:,0],lit[:,1],c='r',label="Literature Values")
#ax.scatter(alex[:,0],alex[:,1],c='b',label="FRAP Values")
#ax.plot(lit[:,0],m*lit[:,0]+c,'k-')

ax.loglog(lit[:,0],lit[:,1],'ro',label="Literature Values")
ax.loglog(alex[:,0],alex[:,1],'bo',label="FRAP Values")
ax.plot(lit[:,0],(lit[:,0]**m*exp(1)**c),'k-')


ax.set_xlabel('molecular weight')
ax.set_ylabel("diffusivity (\SI{}{\micro\meter^2/s})")

plt.legend(loc=1,prop={'size':6})

plt.draw()
plt.savefig('/home/alex_loc/Documents/Talks/PhDSymposium_2015/figs/theresa_data.eps')
plt.show()
raw_input()



