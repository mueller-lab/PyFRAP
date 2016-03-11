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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load embryo object
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
cwd=os.getcwd()
emb=embryo("2D","frap")
emb=emb.load_embryo("/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/A488-dextran-10kDa_20140610/exp1/results/A488-dextran-10kDa_20140610_exp1.pk")
emb.fn_datafolder=emb.fn_datafolder+"/"


#Loading first img
fn=emb.fn_datafolder+emb.file_list[0]

#Load img
data_img = skiio.imread(fn).astype(emb.data_enc)
data_vals=data_img.real
data_vals=data_vals.astype('float')

from mpl_toolkits.axes_grid.axislines import SubplotZero


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting fish
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

#data_vals=emb.im_reg_ICs
#x_grid=linspace(1,emb.data_res_px,emb.data_res_px)
#y_grid=linspace(1,emb.data_res_px,emb.data_res_px)

#x_grid, y_grid= meshgrid(x_grid,y_grid)

##Creating figure
#fig=plt.figure()
#ax=fig.add_subplot(111)

#cflevels=linspace(data_vals.min(),data_vals.max(),100)

#ax.contourf(x_grid,y_grid,data_vals,levels=cflevels)

#ax.get_xaxis().set_visible(False)
#ax.get_yaxis().set_visible(False)

##plt.axis('equal')

##ax.set_xlim([0,511])
##ax.set_ylim([0,511])

#plt.autoscale(enable=True, axis='both', tight=True)

#plt.draw()
#plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/ideal_inhom.png',dpi=200)

#plt.show()
#raw_input("Done")


fig_width_pt = 307.28/2
inches_per_pt = 1.0/72.27
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
	'axes.labelsize': 7,
	'lines.linewidth': 0.4, 
	'lines.markersize'  : 3,            # markersize, in points
          'text.fontsize': 7,
          'legend.fontsize': 7,
         'xtick.labelsize': 7,
          'ytick.labelsize': 7,
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

#Timeseries plot  
fig=plt.figure()
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)

fig.subplots_adjust(bottom=0.25)

for direction in ["right", "top"]:
	ax.axis[direction].set_visible(False)
for direction in ["bottom", "left"]:	
	ax.axis[direction].set_axisline_style("-|>")

#Data	

emb.tvec_data=emb.tvec_data[::20]
emb.out_av_data_pinned_d=emb.out_av_data_pinned_d[::20]
emb.squ_av_data_pinned_d=emb.squ_av_data_pinned_d[::20]

ax.plot(emb.tvec_data,emb.squ_av_data_pinned_d,'r*',label="squ data")
ax.plot(emb.tvec_data,emb.out_av_data_pinned_d,'b*',label="out data")
#ax.plot(emb.tvec_data,emb.slice_av_data_d,'g-')

ax.plot(emb.fits[0].tvec_fit, emb.fits[0].squ_av_fitted_d,'r-',label="squ fitted") 
ax.plot(emb.fits[0].tvec_fit, emb.fits[0].out_av_fitted_d,'b-',label="out fitted") 
#ax.plot(emb.fits[0].tvec_fit, emb.fits[0].slice_av_fitted_d,'g--',label="slice_fitted") 
		

ax.set_xlabel("time (s)")
ax.set_ylabel("intensity (AU)")

lg=plt.legend(bbox_to_anchor=None, loc=4, borderaxespad=0.)
lg.draw_frame(False)

ax.set_yticks([0.0,0.4,0.8,1.2])

plt.autoscale(enable=True, axis='both', tight=True)
plt.draw()
plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/ideal_plot.eps')



