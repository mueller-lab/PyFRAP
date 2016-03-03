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


load_emb=1
analyze=0
simulate=0
pin=0
fit=0

#~~~~~~~~~~~~~~~~~~~~~~~~
#Creating ideal embryo object
#~~~~~~~~~~~~~~~~~~~~~~~~

emb_ideal=embryo("ideal","frap")

#Filepath
cwd=os.getcwd()
emb_ideal.fn_datafolder=cwd+"/testdata/"
emb_ideal.file_list=get_sorted_folder_list(emb_ideal.fn_datafolder,emb_ideal.data_ft)

#For perfomance purpose, take only every third image
emb_ideal.file_list=emb_ideal.file_list[::3]

#Dataset Settings
emb_ideal.nframes=shape(emb_ideal.file_list)[0]
emb_ideal.framerate=30
emb_ideal.tstart=0
emb_ideal.tend=emb_ideal.framerate*(emb_ideal.nframes-1)
emb_ideal.tvec_data=linspace(emb_ideal.tstart,emb_ideal.tend,emb_ideal.nframes)
emb_ideal.steps_data=emb_ideal.nframes
emb_ideal.dt_data=emb_ideal.framerate

#Simulation Settings
emb_ideal.steps_sim=500
emb_ideal.tvec_sim=linspace(emb_ideal.tstart,emb_ideal.tend,emb_ideal.steps_sim)

emb_ideal.fn_mesh=cwd+"/testdata/meshfile.msh"
emb_ideal.fn_geo=cwd+"/testdata/meshfile.geo"
emb_ideal.usemesh=1
emb_ideal.apply_data=0


#~~~~~~~~~~~~~~~~~~~~~~~~
#Creating embryo object
#~~~~~~~~~~~~~~~~~~~~~~~~

emb=embryo("Fish","frap")

#Filepath
cwd=os.getcwd()
emb.fn_datafolder=cwd+"/testdata/"
emb.file_list=get_sorted_folder_list(emb.fn_datafolder,emb.data_ft)

#For perfomance purpose, take only every third image
emb.file_list=emb.file_list[::3]

#Dataset Settings
emb.nframes=shape(emb.file_list)[0]
emb.framerate=30
emb.tstart=0
emb.tend=emb.framerate*(emb.nframes-1)
emb.tvec_data=linspace(emb.tstart,emb.tend,emb.nframes)
emb.steps_data=emb.nframes
emb.dt_data=emb.framerate

emb.volSize_px=5

#Simulation Settings
emb.steps_sim=500
emb.tvec_sim=linspace(emb.tstart,emb.tend,emb.steps_sim)

emb.usemesh=1
emb.fn_mesh=cwd+"/testdata/meshfile.msh"
emb.fn_geo=cwd+"/testdata/meshfile.geo"

if load_emb==1:
	emb_ideal=emb_ideal.load_embryo(cwd+"/testdata/"+"testembryo_ideal.pk")
	emb=emb.load_embryo(cwd+"/testdata/"+"testembryo.pk")
	
if analyze==1:
	#Image Analysis
	emb_ideal = analyze_dataset(emb_ideal)
	emb_ideal.save_embryo(cwd+"/testdata/"+"testembryo_ideal.pk")
	#emb = analyze_dataset(emb)
	#emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")	
	
if simulate==1:
	
	emb_ideal.squ_av_data_d=emb.squ_av_d
	emb_ideal.out_av_data_d=emb.out_av_d
	emb_ideal.slice_av_data_d=emb.slice_av_d
	emb_ideal.tvec_data=emb.tvec_sim
	emb_ideal.steps_data=emb.steps_sim
	emb_ideal.save_embryo(cwd+"/testdata/"+"testembryo_ideal.pk")
	
	#Simulation	
	emb_ideal.plot_cont=0
	emb_ideal.slice_width_px=2
	emb_ideal=diff_reg(emb_ideal)
	emb_ideal.save_embryo(cwd+"/testdata/"+"testembryo_ideal.pk")
	
	#emb.D=10/(emb.conv_fact**2)
	#emb.slice_width_px=1
	#emb.plot_cont=0
	#emb=diff_reg(emb)
	#emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")	
	
if pin==1:
	#emb=pin_conc(emb)
	#emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")
	emb_ideal=pin_conc(emb_ideal)
	emb_ideal.save_embryo(cwd+"/testdata/"+"testembryo_ideal.pk")
	
if fit==1:
	##Fitting
	#emb.fits=[]

	##Fit equ, pinned
	#emb.add_fit(0,"fit equ pinned","default")

	#emb.fits[0].opt_meth="Nelder-Mead"
	#emb.fits[0].fit_prod=1
	#emb.fits[0].fit_degr=0
	#emb.fits[0].fit_pinned=1
	#emb.fits[0].fit_slice=1
	#emb.fits[0].fit_squ=1
	#emb.fits[0].fit_out=0
	#emb.fits[0].UB_D=400
	#emb.fits[0].x0=[10,0,0]
	#emb.fits[0].equ_on=1
	
	
	emb_ideal.fits=[]
	
	#Fit equ, pinned
	emb_ideal.add_fit(0,"fit equ pinned","default")

	emb_ideal.fits[0].opt_meth="Nelder-Mead"
	emb_ideal.fits[0].fit_prod=0
	emb_ideal.fits[0].fit_degr=0
	emb_ideal.fits[0].fit_pinned=1
	emb_ideal.fits[0].fit_slice=1
	emb_ideal.fits[0].fit_squ=1
	emb_ideal.fits[0].fit_out=0
	emb_ideal.fits[0].UB_D=400
	emb_ideal.fits[0].x0=[30,0,0]
	emb_ideal.fits[0].equ_on=1
	
	emb_ideal=dim_fitting(emb_ideal,0)

	
	#emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")
	emb_ideal.save_embryo(cwd+"/testdata/"+"testembryo_ideal.pk")
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#emb.fits[0].plot_fit_pinned()
#emb_ideal.fits[0].plot_fit_pinned()

#emb.fits[0].print_results()
#emb_ideal.fits[0].print_results()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Outputing figures for presentation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Some Settings
from mpl_toolkits.axes_grid.axislines import SubplotZero

fig_width_pt = 307.28/3.2
inches_per_pt = 1.0/72.27
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_height = fig_width*golden_mean      # height in inches
fig_height = 0.87592585241948118
fig_size =  [fig_width,fig_height]


params = {'backend': 'ps',
	'axes.labelsize': 7,
	'lines.linewidth': 0.75, 
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

fig.subplots_adjust(bottom=0.2)

for direction in ["right", "top"]:
	ax.axis[direction].set_visible(False)
for direction in ["bottom", "left"]:	
	ax.axis[direction].set_axisline_style("-|>")

#Data	
#ax.plot(emb_ideal.tvec_data,emb_ideal.squ_av_data_d,'b-')
#ax.plot(emb_ideal.tvec_data,emb_ideal.out_av_data_d,'r-')
#ax.plot(emb_ideal.tvec_data,emb_ideal.slice_av_data_d,'g-')

ax.plot(emb_ideal.tvec_sim,emb_ideal.squ_av_d,'b--',label="squ\_ideal")
#ax.plot(emb_ideal.tvec_sim,emb_ideal.out_av_d,'r--',label="out\_ideal")
ax.plot(emb_ideal.tvec_sim,emb_ideal.slice_av_d,'g--',label="slice\_ideal")

ax.plot(emb.tvec_sim,emb.squ_av_d,'b-.',label="squ\_interpolation")
#ax.plot(emb.tvec_sim,emb.out_av_d,'r-.',label="out\_3D")
ax.plot(emb.tvec_sim,emb.slice_av_d,'g-.',label="slice\_interpolation")

ax.set_xlabel("time (s)")
ax.set_ylabel("intensity (AU)")

lg=plt.legend(bbox_to_anchor=None, loc=4, borderaxespad=0.)
lg.draw_frame(False)

plt.draw()
#plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/compare_ideal.eps')
#plt.show()

fig=plt.figure()
ax=fig.add_subplot(111)

fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(left=0.35)

		
#Ds=[10,emb_ideal.fits[0].D_opt_mu]
Ds=[10,6.3]
prods=[emb.fits[0].prod_opt_scaled,emb_ideal.fits[0].prod_opt_scaled]
names=["inhomogen.","sharp"]

N=2
#ind= arange(N)
ind=asarray([0.5,1.])
width=0.2

bar_D=ax.bar(ind+width, Ds, width, color='0.75')
#ax2=ax.twinx()

#bar_prod=ax2.bar(ind, prods, width, color='c')

#lg=plt.legend([bar_D,bar_prod],["diffusivity","production"],loc=9,bbox_to_anchor=(0.41, 1))
#lg=plt.legend([bar_D],["diffusivity"],loc=9,bbox_to_anchor=(0.5, 1.2),columnspacing=2)
#lg.draw_frame(False)

ax.set_xticks(ind+1.5*width)
ax.set_xticklabels(names)

ax.set_yticks([0,4,8,12,16])
ax.set_xlim([0.6,1.5])
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
#ax.set_ylabel("diffusivity ($\mu\mathrm{m}^2/\mathrm{s}$)")
ax.set_ylabel("diffusivity (\SI{}{\micro\meter^2/s})")
#ax2.set_ylabel("production ($10^{-4}\, 1/s$)")

plt.draw()
plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/bar_ideal.eps')
#plt.show()



#ax.spines['top'].set_color('none')
#ax.xaxis.set_ticks_position('bottom')

#ax.spines['right'].set_color('none')
#ax.yaxis.set_ticks_position('left')





