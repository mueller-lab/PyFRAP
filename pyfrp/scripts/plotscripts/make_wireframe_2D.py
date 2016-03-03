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
emb=emb.load_embryo(cwd+"/testdata/"+"testembryo.pk")

#Loading first img
fn=emb.fn_datafolder+emb.file_list[0]

#Load img
data_img = skiio.imread(fn).astype(emb.data_enc)
data_vals=data_img.real
data_vals=data_vals.astype('float')

down_sampl=4
data_vals=data_vals[0:emb.data_res_px-1:down_sampl,0:emb.data_res_px-1:down_sampl]

res_new=shape(data_vals)[0]

x_grid=linspace(1,emb.data_res_px,res_new)
y_grid=linspace(1,emb.data_res_px,res_new)

x_grid, y_grid= meshgrid(x_grid,y_grid)

cflevels=linspace(data_vals.min()-5,data_vals.max()+5,20)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating wireframe for fish
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
res_wire=20
		
#Define Grid vectors
v_outer=math.acos(((emb.fish_inradius_px**2-emb.fish_outradius_px**2+emb.fish_dist_px**2)/(2*emb.fish_dist_px)-emb.fish_dist_px)/emb.fish_outradius_px)
v_inner=math.acos(((emb.fish_inradius_px**2-emb.fish_outradius_px**2+emb.fish_dist_px**2)/(2*emb.fish_dist_px))/emb.fish_inradius_px)

#Parametric vectors
u=r_[0:2*pi:50j]
v=r_[v_inner:0:50j]
#v=r_[0:pi:50j]

x_wire_inner=emb.center_embr_px[0]+emb.fish_inradius_px*outer(cos(u),sin(v))
y_wire_inner=emb.center_embr_px[1]+emb.fish_inradius_px*outer(sin(u),sin(v))
z_wire_inner=-(emb.fish_dist_px+emb.fish_outradius_px)+emb.fish_inradius_px*outer(ones(size(u)),cos(v))

#Parametric vectors
u=r_[0:2*pi:50j]
v=r_[0:v_outer:50j]
#v=r_[0:pi:50j]

x_wire_outer=emb.center_embr_px[0]+emb.fish_outradius_px*outer(cos(u),sin(v))
y_wire_outer=emb.center_embr_px[1]+emb.fish_outradius_px*outer(sin(u),sin(v))
z_wire_outer=-emb.fish_outradius_px+emb.fish_outradius_px*outer(ones(size(u)),cos(v))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Output settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

from mpl_toolkits.axes_grid.axislines import SubplotZero

fig_width_pt = 307.28/2
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
         'font.sans-serif': 'Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Helvetica, Avant Garde, sans-serif'}#,
         #'figure.figsize': fig_size}
plt.rcParams.update(params)


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Bitstream Vera Sans', 'Lucida Grande', 'Verdana', 'Geneva']  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting fish
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Creating figure
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')

#Finally plot
ax.plot_wireframe(x_wire_outer,y_wire_outer,z_wire_outer,color='b')
ax.plot_wireframe(x_wire_inner,y_wire_inner,z_wire_inner,color='b')

#Plot img
#ax.contourf(x_grid,y_grid,data_vals,levels=cflevels,offset=emb.slice_height_px[0])
		
ax.set_zlim([z_wire_outer.min(),z_wire_outer.max()])

ax.set_frame_on(False)
ax.set_axis_off()

plt.tight_layout()
		
plt.draw()
#plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Outputting fish
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/fish_wire.eps')

#Clearing axes
fig2=plt.figure()
ax=fig2.add_subplot(111,projection='3d')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Load embryo
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

emb=emb.load_embryo(cwd+"/testdata/"+"testembryo_2D.pk")
res_wire=40

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Creating wireframe for plane
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
#Define Grid vectors
x_wire=linspace(-emb.cylinder_radius_px+emb.center_embr_px[0]+0.0001,emb.cylinder_radius_px+emb.center_embr_px[0]-0.0001,res_wire)
z_wire=asarray([0])
X_wire, Z_wire=meshgrid(x_wire,z_wire)
Y_wire=emb.center_embr_px[1]+sqrt(emb.cylinder_radius_px**2-(X_wire-emb.center_embr_px[0])**2)
Y_wire2=-Y_wire+2*emb.center_embr_px[1]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting plane
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
#Finally plot
ax.plot_wireframe(X_wire,Y_wire,Z_wire)
ax.plot_wireframe(X_wire,Y_wire2,Z_wire)

#ax.contourf(x_grid,y_grid,data_vals,levels=cflevels,offset=0)

ax.set_zlim([-5,5])

#ax.grid(False)
ax.set_frame_on(False)
ax.set_axis_off()

plt.draw()
#plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Outputting fish
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

plt.savefig('/home/alex_loc/Documents/Talks/zebrafish2014/figs/2D_wire.eps')



		
		