#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Misc
import os
import time
from optparse import OptionParser

#Numpy/Scipy
from numpy import *
from scipy.interpolate import griddata

#PyFRAP Modules
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_img_module as pyfrp_img
import pyfrp_fit_module as pyfrp_fit
import pyfrp_gmsh_module as pyfrp_gmsh
import pyfrp_plot_module as pyfrp_plot
import pyfrp_integration_module as pyfrp_int
from embryo import *
from molecule import *

#Matplotlib
import matplotlib.pyplot as plt

#Image Processing
import skimage.filter as skifilt

#===========================================================================================================================================================================
#Options
#===========================================================================================================================================================================

nbins=100
bins_fixed=False
density=False
thresh=10
fill=nan
to_steps=800

#===========================================================================================================================================================================
#Script
#===========================================================================================================================================================================

#Load molecule
mol=molecule("bla")
mol=mol.load_molecule("../Ideal_FRAP/Gary/20150731/Analysis.pk")

#Grab first embryo
emb=mol.embryos[1]

#-----------------------------------------------
#Plot concentration over all nodes
fig=plt.figure()
fig.show()
ax=fig.add_subplot(121)
ax.plot(emb.tvec_sim,emb.all_av_d)
ax=fig.add_subplot(122)

print len(emb.tvec_sim[:to_steps])

ax.plot(emb.tvec_sim[:to_steps],emb.squ_av_d[:to_steps],'b')
ax.plot(emb.tvec_sim[:to_steps],emb.out_av_d[:to_steps],'r')
ax.plot(emb.tvec_sim[:to_steps],emb.slice_av_d[:to_steps],'g')

plt.draw()
print min(emb.all_av_d)/max(emb.all_av_d)
raw_input()



#-----------------------------------------------
#Compare different initial images

#Load images
fn_img=emb.fn_datafolder+emb.file_list[0]
img0 = skiio.imread(fn_img).astype(emb.data_enc)		
img0=img0.real
img0=img0.astype('float')

img_pre = skiio.imread(emb.fn_preimage).astype(emb.data_enc)		
img_pre=img_pre.real
img_pre=img_pre.astype('float')

#Gaussian filter
emb.gaussian_sigma=1.5
img0_gauss=skifilt.gaussian_filter(img0,emb.gaussian_sigma)
img_pre_gauss=skifilt.gaussian_filter(img_pre,emb.gaussian_sigma)

#Compute normed image
offset=1
IC=(img0+offset)/(img_pre+offset)
IC_gauss=(img0_gauss+offset)/(img_pre_gauss+offset)

#Mask everything out that is greater than a thresh in IC image
mask_IC=ones(shape(IC))
mask_IC=where(IC<thresh,mask_IC,fill)
print "Pixels in mask_IC = ",nansum(mask_IC)
mask_IC=IC*mask_IC

#Interpolation of mask to fill empty pixels
x=[]
y=[]
z=[]
xint=[]
yint=[]
zint=[]

for i in range(shape(mask_IC)[0]):
	for j in range(shape(mask_IC)[1]):	
		if not isnan(mask_IC[i,j]):
			x.append(i)
			y.append(j)
			z.append(mask_IC[i,j])
		else:
			xint.append(i)
			yint.append(j)

x=asarray(x)
y=asarray(y)
z=asarray(z)

xs=asarray([x,y])

grid_x,grid_y=meshgrid(arange(shape(mask_IC)[0]),arange(shape(mask_IC)[1]))

mask_IC_gauss=griddata((x,y),z,(grid_x,grid_y),method='nearest')

print float(mask_IC_gauss.sum())/float(IC.sum())

print "img0=",pyfrp_int.calc_max_grad(img0)
print "img_pre=",pyfrp_int.calc_max_grad(img_pre)
print "img0_gauss=",pyfrp_int.calc_max_grad(img0_gauss)
print "img_pre_gauss=",pyfrp_int.calc_max_grad(img_pre_gauss)
print "IC=",pyfrp_int.calc_max_grad(IC)
print "IC_gauss=",pyfrp_int.calc_max_grad(IC_gauss)
print "mask_IC_gauss=",pyfrp_int.calc_max_grad(mask_IC_gauss)

raw_input()
	
#Computing histogram
if bins_fixed:
	bins_img=linspace(nanmin(img0),nanmax(img0),256+1)
	bins_gauss=linspace(nanmin(img0_gauss),nanmax(img0_gauss),256+1)
	bins_IC=linspace(nanmin(IC),nanmax(IC),256+1)
	bins_mask=linspace(nanmin(mask_IC),nanmax(mask_IC),256+1)

	hist0,bin_edges_0=histogram(img0,bins_img,density=density)
	hist_pre,bin_edges_pre=histogram(img_pre,bins_img,density=density)

	hist_gauss,bin_edges_gauss=histogram(img0_gauss,bins_gauss,density=density)
	hist_pre_gauss,bin_edges_pre_gauss=histogram(img_pre_gauss,bins_gauss,density=density)

	hist_IC,bin_edges_IC=histogram(IC,bins_IC,density=density)
	hist_IC_gauss,bin_edges_IC_gauss=histogram(IC_gauss,bins_IC,density=density)

	bins_img=linspace(nanmin(img0),nanmax(img0),256)
	bins_gauss=linspace(nanmin(img0_gauss),nanmax(img0_gauss),256)
	bins_IC=linspace(nanmin(IC),nanmax(IC),256)
	bins_mask=linspace(nanmin(mask_IC),nanmax(mask_IC),256)
	
	bins_pre=bins_img.copy()
	bins_pre_gauss=bins_gauss.copy()
	bins_IC_gauss=bins_IC.copy()
	
else:
	hist0,bin_edges_0=histogram(img0,bins=nbins,density=density)
	hist_pre,bin_edges_pre=histogram(img_pre,bins=nbins,density=density)

	hist_gauss,bin_edges_gauss=histogram(img0_gauss,bins=nbins,density=density)
	hist_pre_gauss,bin_edges_pre_gauss=histogram(img_pre_gauss,bins=nbins,density=density)
	
	hist_IC,bin_edges_IC=histogram(IC,bins=nbins,density=density)
	hist_IC_gauss,bin_edges_IC_gauss=histogram(IC_gauss,bins=nbins,density=density)
	
	if isnan(fill):
		bins_mask=linspace(-0.1,thresh,nbins+1)
		hist_mask,bin_edges_mask=histogram(mask_IC,bins_mask,density=density)
	else:		
		hist_mask,bin_edges_mask=histogram(mask_IC,bins=nbins,density=density)
	
	hist_mask_IC_gauss,bin_edges_mask_IC_gauss=histogram(mask_IC_gauss,bins=nbins,density=density)
	
	bins_img=bin_edges_0[:nbins]
	bins_pre=bin_edges_pre[:nbins]
	bins_gauss=bin_edges_gauss[:nbins]
	bins_pre_gauss=bin_edges_pre_gauss[:nbins]
	bins_IC=bin_edges_IC[:nbins]
	bins_IC_gauss=bin_edges_IC_gauss[:nbins]
	bins_mask=bin_edges_mask[:nbins]
	bins_mask_gauss=bin_edges_mask_IC_gauss[:nbins]
		
#Plot
fig=plt.figure()
fig.show()

#img0
ax=fig.add_subplot(261)
ax.set_title('img0')
ax.imshow(img0,vmin=0,vmax=65535)
ax=fig.add_subplot(267)
ax.bar(bins_img,hist0,color='b')

#img_pre
ax=fig.add_subplot(262)
ax.set_title('img_pre')
ax.imshow(img_pre,vmin=0,vmax=65535)
ax=fig.add_subplot(268)
ax.bar(bins_pre,hist_pre,color='b')

#img_gauss
ax=fig.add_subplot(263)
ax.set_title('img_gauss')
ax.imshow(img0_gauss)
ax=fig.add_subplot(2,6,9)
ax.bar(bins_gauss,hist_gauss,color='b')

#img_pre_gauss
ax=fig.add_subplot(264)
ax.set_title('img_pre_gauss')
ax.imshow(img_pre_gauss)
ax=fig.add_subplot(2,6,10)
ax.bar(bins_pre_gauss,hist_pre_gauss,color='b')

#IC
ax=fig.add_subplot(265)
ax.set_title('IC')
ax.imshow(IC)
ax=fig.add_subplot(2,6,11)
ax.bar(bins_IC,hist_IC,color='b')

#IC_gauss
ax=fig.add_subplot(266)
ax.set_title('IC_gauss')
ax.imshow(IC_gauss)
ax=fig.add_subplot(2,6,12)
ax.bar(bins_IC_gauss,hist_IC_gauss,color='b')

plt.draw()

#mask_IC
fig=plt.figure()
fig.show()

ax=fig.add_subplot(221)
ax.set_title('mask_IC')
ax.imshow(mask_IC)
ax=fig.add_subplot(223)
w=(max(bins_mask)-min(bins_mask))/(2*len(bins_mask))
ax.bar(bins_mask,hist_mask,width=w,color='b')

ax=fig.add_subplot(222)
ax.set_title('mask_IC_gauss')
ax.imshow(mask_IC_gauss)
ax=fig.add_subplot(224)
w=(max(bins_mask_gauss)-min(bins_mask_gauss))/(2*len(bins_mask_gauss))
ax.bar(bins_mask_gauss,hist_mask_IC_gauss,width=w,color='b')

plt.draw()
raw_input()


#-----------------------------------------------
#Plot concentration over all nodes
fig=plt.figure()
fig.show()
ax=fig.add_subplot(111)
ax.plot(emb.tvec_sim,emb.all_av_d)
plt.draw()
raw_input()
