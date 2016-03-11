#! /usr/bin/python

from fipy import *
import math 
import os
from numpy import *
from pyfrp_sim_module import *
from pyfrp_misc_module import *
from pyfrp_img_module import *
from pyfrp_fit_module import *
from pyfrp_stats_module import *
from embryo import *
import matplotlib

res_fcs=[141,107,80,64]
res_frap=[87.8216,75.2855,52.3918,37.26145]

res_names=["3kDa water", "3kDa agerose","10kDa water", "10kDa agerose"]

N=shape(res_names)[0]
ind= arange(N)
width=0.3


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

fig=plt.figure()
fig.show()
fig.subplots_adjust(bottom=0.3)
ax=fig.add_subplot(111)

ax.bar(ind-width, res_fcs, width, color='r',label='FCS')
ax.bar(ind, res_frap, width, color='b',label='FRAP')
					
ax.legend(loc=1, borderaxespad=0.)	
		
ax.set_xticks(ind)
ax.set_xticklabels(res_names,rotation=90)

plt.draw()
raw_input()
