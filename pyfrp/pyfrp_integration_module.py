#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Integration module for PyFRAP toolbox, including following functions:

#(1)  avg_conc: Uses arithmetic averaging to compute concentration profiles in selected regions
#(2)  int_conc: Uses integration methods to compute concentration profiles in selected regions
#(3)  calc_conc_int: Calls proper functions for area calculation and integration
#(4)  trapz_int_all: Arranges data for trapezodial integration and offers debugging options, then passes data to either very_fast_trapz_app_int, fast_trapz_app_int or trapz_app_int
#(5)  very_fast_trapz_app_int: Spilts up data in subsections, then uses trapezodial rule for integration
#(6)  fast_trapz_app_int: Uses trapezodial rule for integration without any index mapping
#(7)  trapz_app_int: Uses trapezodial rule for integration with index mapping
#(8)  idx_trapz: Finds index mapping used for trapz_app_int
#(9)  eval_int_trapz: Evaluates integrals using results of idx_trapz and returns it to trapz_app_int
#(10) calc_tet_vol: Calculates tetrahedron volume
#(11) calc_tet_sidelengths: Calculates tetrahedron sidelengths
#(12) calc_slice_vol: Calculates slice volume
#(13) calc_slice_area: Calculates slice area
#(14) calc_slice_img_area: Calculates intersection of slice and image area
#(15) corner_integrand: Integrand needed to compute corner area in calc_slice_img_area
#(16) calc_bleached_vol: Calculates bleached region volume (square assumed)
#(17) calc_bleached_area: Calculates bleached region area (box assumed)

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#from fipy import *
#from numpy import *
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#import matplotlib.collections as collections
#import matplotlib.pyplot as plt
#import matplotlib.patches as ptc
#import matplotlib.tri as tri
#import pylab as plab
#import time
#import sys
#from matplotlib.path import Path
from numpy import linalg as LA
#import scipy.interpolate as interp
#import scipy.integrate as integr
#import scipy.ndimage.interpolation as ndi

#PyFRAP Modules
#import pyfrp_plot_module as pyfrp_plot
#import pyfrp_debug_module as pyfrp_debug


#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Averages over slice and computes concentration profiles for regions

#def avg_conc(embryo,phi,cvs,squ_ind,out_ind,squ_small_ind,out_small_ind,slice_ind,outer_ind,inner_ind,squ_pocket_ind,out_pocket_ind,all_ind,squ_av,out_av,squ_small_av,out_small_av,slice_av,outer_av,inner_av,squ_pocket_av,out_pocket_av,all_av):
	
	##Averaging for large square
	##print phi.value[squ_ind].min(),phi.value[squ_ind].max(), mean(phi.value[squ_ind])
	##raw_input()
	
	
	#squ_av.append(get_avg(phi,cvs,squ_ind))
	#out_av.append(get_avg(phi,cvs,out_ind))
	
	##Averaging for small square
	#if embryo.avg_small==1:
		#squ_small_av.append(get_avg(phi,cvs,squ_small_ind))
		#out_small_av.append(get_avg(phi,cvs,out_small_ind))

	##Averaging for slice
	#slice_av.append(get_avg(phi,cvs,slice_ind))	
	
	##Averaging outer
	#if embryo.avg_outer==1:
		#outer_av.append(get_avg(phi,cvs,outer_ind))
	##Averging inner
	#if embryo.avg_inner==1:
		#inner_av.append(get_avg(phi,cvs,inner_ind))
	##Avering all
	#if embryo.avg_all==1:
		#all_av.append(get_avg(phi,cvs,all_ind))
	
	##Avering pockets
	#if embryo.avg_pocket==1:
		#squ_pocket_av.append(get_avg(phi,cvs,squ_pocket_ind))
		#out_pocket_av.append(get_avg(phi,cvs,out_pocket_ind))
			
	#return squ_av,out_av,squ_small_av,out_small_av,slice_av,outer_av,inner_av,squ_pocket_av,out_pocket_av,all_av

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns 

def getAvgConc(val,cvs,ind):
	if len(ind)>0:
		return sum(val.value[ind]*cvs[ind])/sum(cvs[ind])
	else:
		return 0.

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates tetrahedron volume for 4 given points

#def calc_tet_vol(point0,point1,point2,point3):

	##Taking point0 as base point, calculating vectors
	#vec1=point1-point0
	#vec2=point2-point0
	#vec3=point3-point0
	
	##Calculating volume: V=1/6*| (a x b) * c |
	#vol=1./6.*abs(dot(cross(vec1,vec2),vec3))
	
	#return vol

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculates sidelengths of tetrahedron given by 4 points

def calcTetSidelengths(point0,point1,point2,point3):

	#Taking point0 as base point, calculating vectors
	vec1=point1-point0
	vec2=point2-point0
	vec3=point3-point0
	
	norms=[LA.norm(vec1),LA.norm(vec2), LA.norm(vec3)]
	
	return norms

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculates analytical slice_volume and approximation through cylinder

#def calc_slice_vol(radius,geometry,slice_height,slice_width):
	
	
	#if geometry=="Cylinder":
		##Case Cylinder: Slice volume is a simple cylinder
		#vol=radius**2*pi*2.*slice_width
		
	#elif geometry=="Fish":
		##Case Fish: Slice vol is difference of Ball segments
		##Formala for ball segment volume is V= 1/3 * h^2*pi *(3h-r)
		
		##Larger segment
		#h=abs(slice_height-slice_width)
		#vol_seg_large=h**2*pi/3. * (3.*radius - h)
		
		##Smaller segment
		#h=abs(slice_height+slice_width)
		#vol_seg_small=h**2*pi/3. * (3.*radius - h)
		
		##Difference
		#vol=vol_seg_large-vol_seg_small
	
	#else: 
		#print "ERROR: calc_slice_vol hasn't been implemented for geometry=", geometry
	
	##Calculate approximate slice volume through approximation through cylinder
	#r_at_slice_h=sqrt((radius**2-(radius-abs(slice_height))**2))
	#vol_approx=r_at_slice_h**2*pi*2.*slice_width
	
	#return vol, vol_approx

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates analytical slice_area

#def calc_slice_area(radius):
	
	#area=pi*radius**2

	#return area

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates analytical slice_img_area, that is, the interception of slice and the original image data

#def calc_slice_img_area(radius,center,res):
	
	##Calculate interceptions between circle and img
	#x_inter1=0
	#y_inter2=0
	#x_inter2=-sqrt(radius**2-(y_inter2-center[1])**2)+center[0]
	
	
	##Calculate area of corners via integration
	#area_corner=integr.quad(corner_integrand,x_inter1,x_inter2,args=(radius,center))
	#area_corner=area_corner[0]
	
	##Calculate img area
	#area_img=res**2
	
	##area of image without corners
	#area=area_img-4*area_corner
	
	#return area

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Integrand of corner for calc_slice_img_area

#def corner_integrand(x,radius,center):
	
	#return -sqrt(radius**2-(x-center[0])**2)+center[1]
	
	
##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates analytical bleached volume

#def calc_bleached_vol(sidelength,slice_width):
	
	#vol=slice_width*2.*sidelength**2
	
	#return vol

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates analytical bleached area

#def calc_bleached_area(sidelength):
	
	#area=sidelength**2
	
	#return area

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates maximum difference in image

#def calc_max_grad(img):
	
	#dx=diff(img,axis=0)
	#dy=diff(img,axis=1)
	#d=[dx.max(),dy.max()]
	
	#return max(d)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Integrates over slice and computes concentration profiles for regions

#def int_conc(embryo,phi,squ_av,out_av,squ_small_av,out_small_av,slice_av,L_small,offset_small,int_masks=None,int_masks_all=None):

	##Nodes in slice
	#x_dat=x[slice_ind]
	#y_dat=y[slice_ind]
		
	##phi values in slice 
	#phi_dat=phi.value[slice_ind]
	
	##If no integration masks are given, make them
	#if int_masks==None:
		#int_masks=[]
	#if int_masks_small==None:
		#int_masks_small=[]
			
	##Calculating concentrations in bleached and outside for big square
	#c_ar_slice,c_ar_squ,c_ar_out,int_masks=calc_conc_int(x_dat,y_dat,phi_dat,embryo.int_steps,embryo.radius_embr_px,embryo.center_embr_px,embryo.side_length_bleached_px,embryo.offset_bleached_px,embryo.integration_method,int_masks,embryo.add_rim_sim,res_img,embryo.im_reg_ICs,embryo.img_in_domain,0)
	
	#if embryo.avg_small==1:
		##Calculating concentrations in bleached and outside for small square
		#c_ar_slice,c_ar_squ_small,c_ar_out_small,int_masks_small=calc_conc_int(x_dat,y_dat,phi_dat,embryo.int_steps,embryo.radius_embr_px,embryo.center_embr_px,L_small,offset_small,embryo.integration_method,int_masks_small,embryo.add_rim_sim,res_img,embryo.im_reg_ICs,embryo.img_in_domain,0)
		#out_small_av.append(c_ar_out_small)
		#squ_small_av.append(c_ar_squ_small)
		
	##Appending integrated concentrations
	#out_av.append(c_ar_out)
	#squ_av.append(c_ar_squ)	
	#slice_av.append(c_ar_slice)
	
	#return squ_av,out_av,slice_av, int_masks,int_masks_small

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculates concentrations in slice, squ and out via trapezoidal integration and normation over area

#def calc_conc_int(x_dat,y_dat,phi_dat,int_steps,radius_data,center,sidelength_data,offset_data,fast_opt,masks,add_rim_sim,res,im_reg_ICs,img_in_domain,debug_opt):
	
	#add_corners=add_rim_sim
	
	##Trapezodial Integration
	#V_slice,V_sq,V_out,masks = trapz_int_all(x_dat,y_dat,phi_dat,int_steps,radius_data,center,sidelength_data,offset_data,fast_opt,masks,im_reg_ICs,add_corners,debug_opt)
	
	##Area
	#if img_in_domain==0:
	
		#if add_rim_sim==0:
			#area_slice_an=calc_slice_img_area(radius_data,center,res)
			
		#elif add_rim_sim==1:
			#area_slice_an=calc_slice_area(radius_data)
	
	#elif img_in_domain==1:
		#if add_rim_sim==0:
			#area_slice_an=calc_bleached_area(res)
			
		#elif add_rim_sim==1:
			#area_slice_an=calc_slice_area(radius_data)
	
		
	#area_sq_an=calc_bleached_area(sidelength_data)
	#area_out_an=area_slice_an-area_sq_an
	
	##Normation over area
	#c_ar_slice=V_slice/area_slice_an
	#c_ar_sq=V_sq/area_sq_an
	#c_ar_out=V_out/area_out_an
	
	#return c_ar_slice, c_ar_sq, c_ar_out,masks

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Trapezoidal integration of concentration in slice and bleached area

#def trapz_int_all(x_dat,y_dat,phi_dat,int_steps,radius_data,center,sidelength_data,offset_data,fast_opt,masks,im_reg_ICs,add_corners,debug_opt):
	
	##-------------------------------
	##Steps:
	##(1) interpolate phi_dat over grid (int_steps x intsteps) using griddata
	##(2) Calculate integral (volume) of every grid element using trapz_app_int
	##(3) Find indices of grid elements that are in slice/square respectively using idx_trapz
	##(4) Compute integral in slice and square using eval_int_trapz
	##(5) If debug_opt==1, do some plotting
	##(6) If fast_opt==1, use fast_trapz_app_int to save some time, won't work with debug_opt
	##(6) If fast_opt==2, use very_fast_trapz_app_int to save even more time
	
	##Check input:
	#if debug_opt==1 and fast_opt==1: 
		#print "debug_opt=1 and fast_opt=1 in trapz_int_all doesn't work at the same time, will use debug_opt=0 and fast_opt=1"
		#debug_opt=0
	
	
	#inter_meth="rbf"
	#parts=10
	
	##Interpolation
	#start_time_inter=time.clock()
	
	##Grid Vectors
	#x_vec = linspace(min(x_dat),max(x_dat),int_steps)
	#y_vec = linspace(min(y_dat),max(y_dat),int_steps)
	
	##Putting it into 2D array
	#grid_x, grid_y = meshgrid(x_vec, y_vec)
	#grid_x=grid_x.T
	#grid_y=grid_y.T
	
	#if inter_meth=="griddata":
				
		#grid_f=interp.griddata((x_dat,y_dat),phi_dat,(grid_x,grid_y),method="linear")
	
	#elif inter_meth=="rbf":
		
		##Check if int_steps is odd or even
		#if mod(int_steps,parts)==0:
			#pass
		#else: 
			#if debug_opt==1:
				#print "WARNING: int_steps need to be able to divided up in parts, going to set int_steps=int_steps-mod(int_steps,parts)"
			#int_steps=int_steps-mod(int_steps,parts)	
		
		##Divide grid vectors into smaller vectors
		#x_vecs=[]
		#y_vecs=[]
		
		#int_part=int_steps/parts
		
		#for i in range(parts):
			#x_vecs.append(x_vec[i*int_part:(i+1)*int_part])
			#y_vecs.append(y_vec[i*int_part:(i+1)*int_part])
		
		
		
		##Actual interpolation		
		#intfunc=interp.Rbf(x_dat, y_dat, phi_dat, function='linear')
		
		##Empty vectors
		#grid_f=zeros(shape(grid_x))
		
		##For debugging, interpolate over whole grid at once
		#if debug_opt==1 and parts==2:
			#grid_f_part=intfunc(grid_x,grid_y)
			
		##For debugging, do some figures
		#if debug_opt==1 and parts==2:
			
			#fig=plt.figure()
			#ax=fig.add_subplot(241)
			#ax.contourf(grid_x,grid_y,grid_f_part)
			
		##Interpolating over every part
		#for j in range(shape(y_vecs)[0]):
			#for i in range(shape(x_vecs)[0]):
			
				##Creating tempory grid
				#grid_x_temp, grid_y_temp=meshgrid(x_vecs[i], y_vecs[j])
				
				#grid_x_temp=grid_x_temp.T
				#grid_y_temp=grid_y_temp.T
				
				##Evaluating over partial grids
				#grid_f_temp=intfunc(grid_x_temp,grid_y_temp)
				
				#grid_f[i*int_part:(i+1)*int_part][:,j*int_part:(j+1)*int_part]=grid_f_temp
				
				#if debug_opt==1 and parts==2:
					#idx_fig=245+(i*2)+j
					#ax=fig.add_subplot(idx_fig)
					#ax.contourf(grid_x_temp,grid_y_temp,grid_f_temp)
		
		##For debugging, lot difference between interpolation and img
		#if debug_opt==1 and parts==2:
			#diff_part=grid_f-grid_f_part
			#ax=fig.add_subplot(243)
			#ax.contourf(grid_x,grid_y,diff_part)
			#ax=fig.add_subplot(242)
			#ax.contourf(grid_x,grid_y,grid_f)
			#plt.draw()
			#raw_input()
			

						
	#time_inter=time.clock()-start_time_inter
	
	#if fast_opt==1:
		#V_slice, V_sq,V_out=fast_trapz_app_int(grid_x,grid_y,grid_f,radius_data,center,sidelength_data,offset_data)
		#masks=[]
	#elif fast_opt==2:
		#V_slice, V_sq,V_out,masks=very_fast_trapz_app_int(grid_x,grid_y,grid_f,radius_data,center,sidelength_data,offset_data,masks,im_reg_ICs,add_corners,debug_opt)
	#else:	
		#start_time_trapz=time.clock()
		
		##Volume integration
		#V, grid_x_V, grid_y_V=trapz_app_int(grid_x,grid_y,grid_f)
		
		#time_trapz=time.clock()-start_time_trapz
		#start_time_el=time.clock()
		
		##Finding indices of grid in slice and bleached
		#el_slice_x,el_slice_y,el_sq_x,el_sq_y,el_coords_slice,el_coords_sq=idx_trapz(grid_x,grid_y,radius_data,center,sidelength_data,offset_data)
		
		#time_el=time.clock()-start_time_el
		#start_time_eval=time.clock()
		
		##Calculate integral over slice
		#V_slice,V_vec_slice=eval_int_trapz(el_slice_x,el_slice_y,V)
		
		##Calculate integral over bleached
		#V_sq,V_vec_sq=eval_int_trapz(el_sq_x,el_sq_y,V)
		
		#time_eval=time.clock()-start_time_eval
		
		#masks=[]
		
		##Calculating V_out
		#V_out=V_slice-V_sq
		
	#time_trapz_total=time.clock()-start_time_inter	
	#if debug_opt==1:
		
		#print "-------------------------------------------"
		#print "Integration times:"
		#print "-------------------------------------------"
		
		#print "Interpolation:", time_inter
		#if fast_opt==0:
			#print "Trapezodial:", time_trapz
			#print "Elements:", time_el
			#print "Evaluation:", time_eval
		#print "Total:", time_trapz_total
	
	
	#if debug_opt==1:
		
		
		##Plotting
		#fig_interp=plt.figure()
		#ax_interp=fig_interp.add_subplot(121,projection='3d')
		
		#ax_interp.plot_surface(grid_x,grid_y,grid_f,cmap='jet')
		#surf_interp = ptc.Rectangle((0, 0), 1, 1, fc="b")
		
		#ax_interp.plot(x_dat,y_dat,phi_dat,'or')
		#scat_interp=ptc.Circle((0,0),1,fc="r")
		
		#ax_integr=fig_interp.add_subplot(122,projection='3d')
		#if fast_opt==0:
			
			#ax_integr.plot_surface(grid_x,grid_y,V,cmap='jet')
			
		##Titles
		#ax_interp.set_title("Interpolation")
		#ax_integr.set_title("Integral over grid")
		
		##Labels
		#ax_interp.set_xlabel("x")
		#ax_integr.set_xlabel("x")
		
		#ax_interp.set_ylabel("y")
		#ax_integr.set_ylabel("y")
		
		#ax_interp.set_zlabel("phi")
		#ax_integr.set_zlabel("int(phi)")
		
		##Legends
		#ax_interp.legend([surf_interp,scat_interp],['phi_int',"phi_dat"])
		
		#plt.show()
		
		
		#raw_input()
	
	#return V_slice,V_sq,V_out, masks

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Super fast version of trapz_app_int

#def very_fast_trapz_app_int(grid_x,grid_y,grid_f,radius,center,sidelength,offset,masks,im_reg_ICs,add_corners,debug_opt):
	
	##Increments for REGULAR grid
	#dx=grid_x[1,0]-grid_x[0,0]
	#dy=grid_y[0,1]-grid_y[0,0]
	
	#res=shape(im_reg_ICs)[0]
	
	##Creating masks for out, squ and slice
	#if shape(masks)[0]==0:
		
		#mask_slice=zeros(shape(grid_f))
		#mask_out=zeros(shape(grid_f))
		#mask_squ=zeros(shape(grid_f))
		
		##Go through all interpolation coordinates
		#for i in range(shape(grid_x)[0]-1):
			#for j in range(shape(grid_y)[1]-1):
				##Check if in slice
				#if sqrt((grid_x[i,j]-center[0])**2+(grid_y[i,j]-center[1])**2)<radius:
					
					##Check if img completely in domain, then we actually need to check if there is 
					#if add_corners==1:
						#mask_slice[i,j]=1
						
						##Check if in squ
						#if offset[0]<=grid_x[i,j] and offset[1]<=grid_y[i,j] and grid_x[i,j]<=offset[0]+sidelength and grid_y[i,j]<=offset[1]+sidelength:
							#mask_squ[i,j]=1	
						#else:	
							#mask_out[i,j]=1
							
					#elif add_corners==0:
						
						#if 0<=grid_x[i,j] and 0<=grid_y[i,j] and grid_x[i,j]<=res and grid_y[i,j]<=res:
							#mask_slice[i,j]=1
						
							##Check if in squ
							#if offset[0]<=grid_x[i,j] and offset[1]<=grid_y[i,j] and grid_x[i,j]<=offset[0]+sidelength and grid_y[i,j]<=offset[1]+sidelength:
								#mask_squ[i,j]=1	
							#else:	
								#mask_out[i,j]=1
		
		##Build mask array
		#masks=[mask_slice,mask_squ,mask_out]
	
	##If masks are already given, use them
	#else:
		#mask_slice=masks[0]
		#mask_squ=masks[1]
		#mask_out=masks[2]
	
	##Some debugging
	#if debug_opt==1:
		
		#V_slice=grid_f*mask_slice
		#V_squ=grid_f*mask_squ
		#V_out=grid_f*mask_out
		
		#fig=plt.figure()
		
		#min_img=amin(im_reg_ICs)
		#max_img=amax(im_reg_ICs)
		
		#ax_img=fig.add_subplot(241)
		#if shape(im_reg_ICs)[0]==shape(grid_f)[0]:
			#plt_img=ax_img.contourf(mask_slice*im_reg_ICs,vmin=min_img,vmax=max_img)
			#circle_data1=ptc.Circle(center,radius=radius,fill=False,linewidth=3,color='r')
			#ax_img.add_patch(circle_data1)
		
		#else:
			#plt_img=ax_img.contourf(im_reg_ICs,vmin=min_img,vmax=max_img)
			#circle_data1=ptc.Circle(center,radius=radius,fill=False,linewidth=3,color='r')
			#ax_img.add_patch(circle_data1)
		
		#fig.colorbar(plt_img, ticks=[min_img,(max_img-min_img)/2,max_img], orientation='horizontal')
		
		#ax_int_slice=fig.add_subplot(242)
		#plt_int_slice=ax_int_slice.contourf(grid_x,grid_y,V_slice,vmin=min_img,vmax=max_img)
		#fig.colorbar(plt_int_slice, ticks=[min_img,(max_img-min_img)/2,max_img], orientation='horizontal')
		#circle_data2=ptc.Circle(center,radius=radius,fill=False,linewidth=3,color='r')
		#ax_int_slice.add_patch(circle_data2)
		
		#ax_int_squ=fig.add_subplot(243)
		#plt_int_squ=ax_int_squ.contourf(grid_x,grid_y,V_squ,vmin=min_img,vmax=max_img)
		#fig.colorbar(plt_int_squ, ticks=[min_img,(max_img-min_img)/2,max_img], orientation='horizontal')
		#circle_data3=ptc.Circle(center,radius=radius,fill=False,linewidth=3,color='r')
		#ax_int_squ.add_patch(circle_data3)
		
		#ax_int_out=fig.add_subplot(244)
		#plt_int_out=ax_int_out.contourf(grid_x,grid_y,V_out,vmin=min_img,vmax=max_img)
		#fig.colorbar(plt_int_out, ticks=[min_img,(max_img-min_img)/2,max_img], orientation='horizontal')
		#circle_data4=ptc.Circle(center,radius=radius,fill=False,linewidth=3,color='r')
		#ax_int_out.add_patch(circle_data4)
		
		#ax_mask_slice=fig.add_subplot(246)
		#plt_mask_slice=ax_mask_slice.contourf(grid_x,grid_y,mask_slice,vmin=0,vmax=1)
		#fig.colorbar(plt_mask_slice, ticks=[0,1/2,1], orientation='horizontal')
		
		#ax_mask_squ=fig.add_subplot(247)
		#plt_mask_squ=ax_mask_squ.contourf(grid_x,grid_y,mask_squ,vmin=0,vmax=1)
		#fig.colorbar(plt_mask_squ, ticks=[0,1/2,1], orientation='horizontal')
		
		#ax_mask_out=fig.add_subplot(248)
		#plt_mask_out=ax_mask_out.contourf(grid_x,grid_y,mask_out,vmin=0,vmax=1)
		#fig.colorbar(plt_mask_out, ticks=[0,1/2,1], orientation='horizontal')
		
		
		##Plotting differences between image and interpolation
		#if shape(im_reg_ICs)[0]==shape(grid_f)[0]:
			#V_slice_img=im_reg_ICs*mask_slice
			#diff_int_img=100*abs(V_slice_img-V_slice)/V_slice_img
			
			#ax_diff=fig.add_subplot(245)
			#plt_diff=ax_diff.contourf(grid_x,grid_y,diff_int_img,vmin=0,vmax=100)
			#fig.colorbar(plt_diff, orientation='horizontal')
			
			#print nanmax(diff_int_img)
		
		#plt.draw()
	
		
		##Printing out sums (only works with int_steps==res)
		#if shape(im_reg_ICs)[0]==shape(grid_f)[0]:
			##Sums over img
			#V_slice_img=nansum(im_reg_ICs*mask_slice)
			#V_squ_img=nansum(im_reg_ICs*mask_squ)
			#V_out_img=nansum(im_reg_ICs*mask_out)
			
			##Sums over interpolation
			#V_slice=nansum(V_slice)
			#V_squ=nansum(V_squ)
			#V_out=nansum(V_out)
		
			#print "V_slice = ", V_slice, " V_slice_img = ", V_slice_img, " diff = ", abs(V_slice-V_slice_img), " diff in % = ", 100*abs(V_slice-V_slice_img)/V_slice_img 
			#print "V_squ = ", V_squ, " V_squ_img = ", V_squ_img, " diff = ", abs(V_squ-V_squ_img), " diff in % = ", 100*abs(V_squ-V_squ_img)/V_squ_img 
			#print "V_out = ", V_out, " V_out_img = ", V_out_img, " diff = ", abs(V_out-V_out_img), " diff in % = ", 100*abs(V_out-V_out_img)/V_out_img 
			
		#raw_input()
	
	##Overwrite debugging values
	#V_slice=nansum(dx*dy*grid_f*mask_slice)
	#V_squ=nansum(dx*dy*grid_f*mask_squ)
	#V_out=nansum(dx*dy*grid_f*mask_out)
			
	#return V_slice,V_squ,V_out, masks

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Faster version of trapz_app_int

#def fast_trapz_app_int(grid_x,grid_y,grid_f,radius,center,sidelength,offset):
	
	##Increments for REGULAR grid
	#dx=grid_x[1,0]-grid_x[0,0]
	#dy=grid_y[0,1]-grid_y[0,0]
	
	##Empty volumes
	#V_slice=0
	#V_out=0
	#V_squ=0
	
	#for i in range(shape(grid_x)[0]-1):
		#for j in range(shape(grid_y)[1]-1):
			
			##Midpoint function value
			#f_mid=1./4. * (grid_f[i,j]+grid_f[i+1,j]+grid_f[i,j+1]+grid_f[i+1,j+1])
			
			##Calculating Volume
			#V=f_mid*dx*dy
			
			##Check to which integral V is added to
			#if isnan(V):
				#pass
			#else:	
				#if sqrt((grid_x[i,j]-center[0])**2+(grid_y[i,j]-center[1])**2)<radius:
					#V_slice=V_slice+V
					
					#if offset[0]<=grid_x[i,j] and offset[1]<=grid_y[i,j] and grid_x[i,j]<=offset[0]+sidelength and grid_y[i,j]<=offset[1]+sidelength:
						#V_squ=V_squ+V	
					#else:	
						#V_out=V_out+V
				
	#return V_slice,V_squ,V_out

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Approximate integral over grid by calculating volumes over single cells. 

#def trapz_app_int(grid_x,grid_y,grid_f):
	
	##Empty matrix
	#V=zeros((shape(grid_f)[0]-1,shape(grid_f)[1]-1))
	#grid_x_V=zeros(shape(V))
	#grid_y_V=zeros(shape(V))
	
	#for i in range(shape(grid_x)[0]-1):
		#for j in range(shape(grid_x)[1]-1):
			
			##Increments
			#dx=grid_x[i+1,j]-grid_x[i,j]
			#dy=grid_y[i,j+1]-grid_y[i,j]
			
			##Midpoint function value
			#f_mid=1./4. * (grid_f[i,j]+grid_f[i+1,j]+grid_f[i,j+1]+grid_f[i+1,j+1])
			
			##Calculating Volume
			#V[i,j]=f_mid*dx*dy
			
			##Generating grid for plotting
			#grid_x_V[i,j]=grid_x[i,j]+dx/2.
			#grid_y_V[i,j]=grid_y[i,j]+dy/2.
			
	#return V, grid_x_V, grid_y_V

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Find idx for slice and sq 

#def idx_trapz(grid_x,grid_y,radius,center,sidelength,offset):
	
	##Empty lists
	#el_slice_x=[]
	#el_slice_y=[]
	#el_sq_x=[]
	#el_sq_y=[]
	
	#el_coords_slice=zeros((shape(grid_x)[0]-1,shape(grid_x)[1]-1))
	#el_coords_sq=zeros((shape(grid_x)[0]-1,shape(grid_x)[1]-1))
	
	#for i in range(shape(grid_x)[0]-1):
		#for j in range(shape(grid_x)[1]-1):
			#r=sqrt((grid_x[i,j]-center[0])**2+(grid_y[i,j]-center[1])**2)
			
			#if r<radius:
				#el_slice_x.append(i)
				#el_slice_y.append(j)
				
				#el_coords_slice[i,j]=1
				
			#if offset[0]<=grid_x[i,j] and offset[1]<=grid_y[i,j] and grid_x[i,j]<=offset[0]+sidelength and grid_y[i,j]<=offset[1]+sidelength:
				#el_sq_x.append(i)
				#el_sq_y.append(j)
				
				#el_coords_sq[i,j]=1
				
	#return el_slice_x, el_slice_y, el_sq_x, el_sq_y, el_coords_slice, el_coords_sq

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Evaluate integral on idxs

#def eval_int_trapz(idxs_x,idxs_y,V):
	
	#V_sum=0
	#V_vec=[]
	
	#for i in range(shape(idxs_x)[0]):
			
		#idx_x=idxs_x[i]
		#idx_y=idxs_y[i]
		
		#if isnan(V[idx_x,idx_y]):
			#pass
		#else:
			#V_sum=V_sum+V[idx_x,idx_y]
		
		#V_vec.append(V_sum)
	
	#return V_sum, V_vec	
