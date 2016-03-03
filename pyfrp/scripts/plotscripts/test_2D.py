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
#Creating 2D embryo object
#~~~~~~~~~~~~~~~~~~~~~~~~

emb_2D=embryo("2D","frap")

#Filepath
cwd=os.getcwd()
emb_2D.fn_datafolder=cwd+"/testdata/"
emb_2D.file_list=get_sorted_folder_list(emb_2D.fn_datafolder,emb_2D.data_ft)

#For perfomance purpose, take only every third image
emb_2D.file_list=emb_2D.file_list[::3]

#Dataset Settings
emb_2D.nframes=shape(emb_2D.file_list)[0]
emb_2D.framerate=30
emb_2D.tstart=0
emb_2D.tend=emb_2D.framerate*(emb_2D.nframes-1)
emb_2D.tvec_data=linspace(emb_2D.tstart,emb_2D.tend,emb_2D.nframes)
emb_2D.steps_data=emb_2D.nframes
emb_2D.dt_data=emb_2D.framerate

#Simulation Settings
emb_2D.steps_sim=500
emb_2D.tvec_sim=linspace(emb_2D.tstart,emb_2D.tend,emb_2D.steps_sim)

emb_2D.geometry="Cylinder"
emb_2D.cylinder_radius_mu=emb_2D.radius_embr_mu
emb_2D.cylinder_height_mu=10
emb_2D.cylinder_radius_px=emb_2D.cylinder_radius_mu/emb_2D.conv_fact
emb_2D.cylinder_height_px=10
emb_2D.slice_height_px=[-5.]
emb_2D.fn_mesh=cwd+"/testdata/meshfile_2d.msh"
emb_2D.fn_geo=cwd+"/testdata/meshfile_2d.geo"
emb_2D.usemesh=1


update_parm_in_file(emb_2D.fn_geo,"radius",emb_2D.radius_embr_px)
update_parm_in_file(emb_2D.fn_geo,"height",emb_2D.cylinder_height_px)
update_parm_in_file(emb_2D.fn_geo,"center_x",emb_2D.center_embr_px[0])
update_parm_in_file(emb_2D.fn_geo,"center_y",emb_2D.center_embr_px[1])
update_parm_in_file(emb_2D.fn_geo,"volSize_px",emb_2D.volSize_px)

#os.system("gmsh -3 " + emb_2D.fn_geo)

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



#update_parm_in_file(emb.fn_geo,"radius",emb_2D.radius_embr_px)
#update_parm_in_file(emb.fn_geo,"height",emb_2D.cylinder_height_px)
#update_parm_in_file(emb.fn_geo,"center_x",emb_2D.center_embr_px[0])
#update_parm_in_file(emb.fn_geo,"center_y",emb_2D.center_embr_px[1])
#update_parm_in_file(emb.fn_geo,"volSize",emb.volSize_px)


if load_emb==1:
	emb_2D=emb_2D.load_embryo(cwd+"/testdata/"+"testembryo_2D.pk")
	emb=emb.load_embryo(cwd+"/testdata/"+"testembryo.pk")
	
	emb.usemesh=1
	emb.fn_mesh=cwd+"/testdata/meshfile.msh"
	emb.fn_geo=cwd+"/testdata/meshfile.geo"
	
if analyze==1:
	#Image Analysis
	emb_2D = analyze_dataset(emb_2D)
	emb_2D.save_embryo(cwd+"/testdata/"+"testembryo_2D.pk")
	emb = analyze_dataset(emb)
	emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")
	
	print 
	print "Done with data analysis"
	
if simulate==1:
	#Simulation	
	#emb_2D=diff_reg(emb_2D)
	#emb_2D.save_embryo(cwd+"/testdata/"+"testembryo_2D.pk")
	emb.volSize_px=15
	emb.D=10/(emb.conv_fact**2)
	#emb=diff_reg(emb)
	#emb.save_embryo(cwd+"/testdata/"+"testembryo.pk")
	
	emb_2D.squ_av_data_d=emb.squ_av_d
	emb_2D.out_av_data_d=emb.out_av_d
	emb_2D.slice_av_data_d=emb.slice_av_d
	emb_2D.tvec_data=emb.tvec_sim
	emb_2D.steps_data=emb.steps_sim
	
	
	
	print 
	print "Done with simulation"

if pin==1:
	
	emb_2D=pin_conc(emb_2D)
	emb_2D.save_embryo(cwd+"/testdata/"+"testembryo_2D.pk")
	
if fit==1:
	
	emb_2D.fits=[]
	
	#Fit equ, pinned
	emb_2D.add_fit(0,"fit equ pinned","default")

	emb_2D.fits[0].opt_meth="Nelder-Mead"
	emb_2D.fits[0].fit_prod=0
	emb_2D.fits[0].fit_degr=0
	emb_2D.fits[0].fit_pinned=1
	emb_2D.fits[0].fit_slice=1
	emb_2D.fits[0].fit_out=0
	emb_2D.fits[0].fit_squ=1
	emb_2D.fits[0].UB_D=400
	emb_2D.fits[0].x0=[50,0,0]
	emb_2D.fits[0].equ_on=1
	
	emb_2D=dim_fitting(emb_2D,0)
	emb=dim_fitting(emb,0)

	emb_2D.save_embryo(cwd+"/testdata/"+"testembryo_2D.pk")
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#emb.fits[0].plot_fit_pinned()
#emb_2D.fits[0].plot_fit_unpinned()


#emb_2D.fits[0].print_results()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Outputting figures for presentation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


fig=plt.figure()


#ax=fig.add_subplot(111)
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)

fig.subplots_adjust(bottom=0.2)
#for direction in ["xzero", "yzero"]:
	#ax.axis[direction].set_axisline_style("-|>")
	#ax.axis[direction].set_visible(True)

for direction in ["right", "top"]:
	ax.axis[direction].set_visible(False)
for direction in ["bottom", "left"]:	
	ax.axis[direction].set_axisline_style("-|>")
	
#ax.plot(emb_2D.tvec_data,emb_2D.squ_av_data_d,'b-')
#ax.plot(emb_2D.tvec_data,emb_2D.out_av_data_d,'r-')
#ax.plot(emb_2D.tvec_data,emb_2D.slice_av_data_d,'g-')

ax.plot(emb_2D.tvec_sim,emb_2D.squ_av_d,'b--',label="squ\_2D")
#ax.plot(emb_2D.tvec_sim,emb_2D.out_av_d,'r--',label="out\_2D")
ax.plot(emb_2D.tvec_sim,emb_2D.slice_av_d,'g--',label="slice\_2D")

ax.plot(emb.tvec_sim,emb.squ_av_d,'b-.',label="squ\_3D")
#ax.plot(emb.tvec_sim,emb.out_av_d,'r-.',label="out\_3D")
ax.plot(emb.tvec_sim,emb.slice_av_d,'g-.',label="slice\_3D")

ax.set_xlabel("time (s)")
ax.set_ylabel("intensity (AU)")

lg=plt.legend(bbox_to_anchor=None, loc=4, borderaxespad=0.)

lg.draw_frame(False)


#ax.spines['top'].set_color('none')
#ax.xaxis.set_ticks_position('bottom')

#ax.spines['right'].set_color('none')
#ax.yaxis.set_ticks_position('left')




#ax.axis["top"].set_visible(False)

#for direction in ["xzero", "yzero"]:
	#ax.axis[direction].set_axisline_style("-|>")
	#ax.axis[direction].set_visible(True)

#plt.draw()
#plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/compare_2D.eps')
#plt.show()	

fig=plt.figure()
ax=fig.add_subplot(111)

fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(left=0.35)

		
Ds=[10,emb_2D.fits[0].D_opt_mu]
#prods=[emb.fits[0].prod_opt_scaled,emb_2D.fits[0].prod_opt_scaled]
names=["3D","2D"]

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
plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/bar_2D.eps')
#plt.show()



	