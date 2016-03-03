#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Debugging module for image analysis and simulation for PyFRAP toolbox, including following functions:
#(1)  debug_geometry: Some debugging messages for geometry.
#(2)  debug_mesh: Some debugging messages for mesh properties.
#(3)  debug_slices: Some debugging messages for chosen slices.
#(4)  debug_squ_ind: Some debugging messages for indices found for squ/out regions.
#(5)  debug_averaging: Some debugging messages to compare both averaging methods.
#(6)  debug_IC_interp: Interpolates over all given slice nodes and returns plots/errors for each slice.

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
from numpy import *

#matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Misc
import sys

#PyFRAP
import pyfrp_integration_module as pyfrp_integr

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out geometry properties of embryo object

def debug_geometry(embryo):
	
	print "-------------------------------------------"
	print "Geometry:" ,embryo.geometry
	print "-------------------------------------------"
		
	#Loading embryo.geometry settings
	if embryo.geometry=="Cylinder":
		slice_range=[-embryo.cylinder_height_px,0]
		
		print "embryo.cylinder_height_px=", embryo.cylinder_height_px
		print "embryo.cylinder_radius_px=", embryo.cylinder_radius_px
		
	elif embryo.geometry=="Fish":	
		
		#Calculating additional pts
		z_interc=-((embryo.fish_outradius_px**2-embryo.fish_inradius_px**2)/(2*embryo.fish_dist_px) - embryo.fish_outradius_px + embryo.fish_dist_px/2)
		x_interc=sqrt(embryo.fish_inradius_px**2-((embryo.fish_outradius_px**2-embryo.fish_inradius_px**2)/embryo.fish_dist_px**2))
		x_interc=z_interc
		
		print "embryo.fish_outradius_px=",embryo.fish_outradius_px
		print "embryo.fish_inradius_px=",embryo.fish_inradius_px
		print "embryo.fish_dist_px=", embryo.fish_dist_px
		
		print "Intersection of balls at (r,z)=(",x_interc,",",z_interc,")"
		
		print "z_intersection-embryo.fish_outradius_px (should be =0) =",z_interc-embryo.fish_outradius_px
		print "x_intersection-embryo.fish_outradius_px (should be =0) =", x_interc-embryo.fish_outradius_px
		
		slice_range=[-embryo.fish_outradius_px-embryo.fish_dist_px+embryo.fish_inradius_px,0]
		
	elif embryo.geometry=="Frog":
		print "embryo.frog_radius_px", embryo.frog_radius_px
		slice_range=[-2*embryo.frog_radius_px,0]
		
	print "Slices are allowed in:", slice_range
	
	return slice_range

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out mesh properties

def debug_mesh(embryo):
	
	#Extracting mesh data
	cellvols=embryo.mesh.getCellVolumes()
	ordvert=embryo.mesh._getOrderedCellVertexIDs()

	#Calculating sidelengths of tetrahedrons
	sls_vec=[]
	for i in range(shape(ordvert)[1]):
		curr_vert=ordvert[:,i]
		
		pt1=embryo.mesh.vertexCoords[:,curr_vert[0]]
		pt2=embryo.mesh.vertexCoords[:,curr_vert[1]]
		pt3=embryo.mesh.vertexCoords[:,curr_vert[2]]
		pt4=embryo.mesh.vertexCoords[:,curr_vert[3]]
		
		sl1,sl2,sl3=pyfrp_integr.calc_tet_sidelengths(pt1,pt2,pt3,pt4)
		
		sls_vec.append(sl1)
		sls_vec.append(sl2)
		sls_vec.append(sl3)
		
	sls_avg=mean(sls_vec)
	sls_max=max(sls_vec)
	sls_min=min(sls_vec)

	print "-------------------------------------------"
	print "Mesh:"
	print "-------------------------------------------"
	print "Mesh has ", shape(embryo.mesh.x)[0] , " cells"
	print "Mesh has ", embryo.mesh._numberOfVertices, "vertices"
	print "Mesh has ", embryo.mesh.numberOfFaces, "faces"
	print "min x=", min(embryo.mesh.x), "max x=", max(embryo.mesh.x)
	print "min y=", min(embryo.mesh.y), "max y=", max(embryo.mesh.y)
	print "min z=", min(embryo.mesh.z), "max z=", max(embryo.mesh.z)
	print "Maximum cell volume= ", max(embryo.mesh.getCellVolumes())
	print "Minimum cell volume= ", min(embryo.mesh.getCellVolumes())
		
	print "Maximum cell volume is", max(cellvols), "in cell number=", argmax(cellvols)
	print "Minimum cell volume is", min(cellvols), "in cell number=", argmin(cellvols)
	print "Average cell volume is", mean(cellvols)

	print "Average sidelength of tetrahedron in embryo.mesh:", sls_avg
	print "Maximum sidelength of tetrahedron in embryo.mesh:", sls_max
	print "Minimum sidelength of tetrahedron in embryo.mesh:", sls_min	
	
	return

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out mesh properties
	
def debug_slices(embryo,slice_range):
	
	#Check if there are no or more than one slice
	if shape(embryo.slice_height_px)[0]<1:
		print "There is no slice given, this will cause problems!"
	elif shape(embryo.slice_height_px)[0]>1:
		print "There are more than one slice given, going to only use first slice for calculation of concentration profiles!"
	
	slice_remove=[]
	
	for j in range(shape(embryo.slice_height_px)[0]):
	
		#Check if slice is in the right range
		if min(slice_range)>embryo.slice_height_px[j] or embryo.slice_height_px[j]>max(slice_range):	
			print "WARNING: Slice ", j , "has height z= ", embryo.slice_height_px[j], "while range of diffusion domain is", min(slice_range) , " - ", max(slice_range), "." 
			print "-> you might encounter problems with your concentration profiles. For now, I am going to remove this slice"
			
			slice_remove.append(embryo.slice_height_px[j])
		
	#Removing slices that are out of domain	
	if shape(slice_remove)[0]>0:
		for entr in slice_remove:
			embryo.slice_height_px.remove(entr)	
		
	return embryo

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out debug messages for region indices

def debug_squ_ind(embryo,phi,L,offset,squ_ind,out_ind,inner_ind,outer_ind,out_pocket_ind,squ_pocket_ind,out_pocket_ind_img,out_ind_img):
	
	x,y,z=embryo.mesh.cellCenters
	
	
	print "-------------------------------------------"
	print "Finding square entries with L=",L, " offset=", offset, " slice_height_px=", embryo.slice_height_px, " slice_width_px=", embryo.slice_width_px, " slice_bottom=", embryo.slice_bottom, ":"
	print "-------------------------------------------"
	
	if shape(squ_ind)[0]<1:
		print "WARNING: There is no node in squ!!!"
	else:
		print "Square sidelength=", embryo.side_length_bleached_px
		print "Square reaches from:", embryo.center_embr_px[0]-L/2 ,"<x,y<" , embryo.center_embr_px[0]+L/2 , "using embryo.center_embr_px"
		print "Square reaches from:", embryo.offset_bleached_px[0] ,"<x,y<" , embryo.offset_bleached_px[0]+embryo.side_length_bleached_px , "using offset"
		
		if embryo.slice_bottom==0:
			print "Slice reaches from",embryo.slice_height_px[0]-embryo.slice_width_px, "<z<" , embryo.slice_height_px[0]+embryo.slice_width_px
		elif embryo.slice_bottom==1:
			print "Slice reaches from",embryo.slice_height_px[0], "<z<" , embryo.slice_height_px[0]+2*embryo.slice_width_px	
		
		print shape(embryo.mesh.x), max(squ_ind), min(squ_ind) 
		print "max(x)=", max(x[squ_ind]), " min(x)=", min(x[squ_ind])
		print "max(y)=", max(y[squ_ind]), " min(y)=", min(y[squ_ind])
		print "max(z)=", max(z[squ_ind]), " min(z)=", min(z[squ_ind])
		print "max(phi)=", max(phi.value[squ_ind]), " min(phi)=", min(phi.value[squ_ind])
		print "With embryo.add_rim_sim=1:","#out=",shape(out_ind)[0] , " #squ=", shape(squ_ind)[0], " #slice=", shape(squ_ind)[0]+shape(out_ind)[0]
		print "With embryo.add_rim_sim=0:","#out=",shape(out_ind_img)[0] , " #squ=", shape(squ_ind)[0], " #slice=", shape(squ_ind)[0]+shape(out_ind_img)[0]

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out debug messages for comparison of the two averaging methods
	
def debug_averaging(embryo,squ_av,out_av,slice_av,squ_av_old,out_av_old,slice_av_old):
		
	print "-------------------------------------------"
	print "Averaging test:"
	print "-------------------------------------------"
	
	print "Slice via interpolation over area:", slice_av[0], ",via avg=", slice_av_old[0], ", data:", embryo.slice_av_data_d[0],"avg err in %:", abs(slice_av_old[0]-embryo.slice_av_data_d[0])/embryo.slice_av_data_d[0],"int err in %:", abs(c_ar_slice-embryo.slice_av_data_d[0])/embryo.slice_av_data_d[0]
	print "Squ via interpolation over area:", squ_av[0], ",via avg=", squ_av_old[0], ", data:", embryo.squ_av_data_d[0],"avg err in %:", abs(squ_av_old[0]-embryo.squ_av_data_d[0])/embryo.squ_av_data_d[0],"int err in %:", abs(c_ar_squ-embryo.squ_av_data_d[0])/embryo.squ_av_data_d[0]
	print "Out via interpolation over area:", out_av_old[0], ",via avg=", out_av_old[0], ", data:", embryo.out_av_data_d[0],"avg err in %:", abs(out_av_old[0]-out_av_data_d[0])/out_av_data_d[0],"int err in %:", abs(c_ar_out-out_av_data_d[0])/out_av_data_d[0]
			
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Interpolate over phi to check ICs if ICs are applied properly

def debug_interpolation_appl(x,y,phi,slice_inds,im_reg_ICs,radius_data,center):
	
	#Getting img resolution
	res=shape(im_reg_ICs)[0]
	
	slice_ind=slice_inds[0]
	
	#Getting slice enties of x,y,phi
	phi_slice=phi.value[slice_ind]
	x_slice=x[slice_ind]
	y_slice=y[slice_ind]
	
	#Interpolation
	phi_int=interp.Rbf(x_slice, y_slice, phi_slice, function='linear')
	
	#Creating interpolation grid
	#x_int = linspace(min(x_slice), max(x_slice), 200)
	#y_int = linspace(min(y_slice), max(y_slice), 200)
	x_int = linspace(0, 512, 200)
	y_int = linspace(0, 512, 200)
	x_int_grid,y_int_grid=meshgrid(x_int,y_int)
	
	#Mapping interpolation to interpolation grid
	phi_slice_int=zeros((len(x_int),len(y_int)))
	
	phi_int_max=0.
	phi_int_min=inf		
	
	print
	print "Interpolating over first slice:"
	
	for i in range(len(x_int)):
		for j in range(len(y_int)):
			
			phi_slice_int[i,j]=phi_int(x_int[i],y_int[j])
			
			if phi_int_min>phi_slice_int[i,j]:
				phi_int_min=phi_slice_int[i,j]
			elif phi_int_max<phi_slice_int[i,j]:
				phi_int_max=phi_slice_int[i,j]
		
		#Print Progress
		curr_perc=int(100*i/float(shape(x_int)[0]))
		sys.stdout.write("\r%d%%" %curr_perc)  
		sys.stdout.flush()
	
	
	#Interpolation on same grid as img for difference calculation
	x_int_diff=linspace(1, 512, 512)
	y_int_diff=linspace(1, 512, 512)
	
	x_int_diff_grid,y_int_diff_grid=meshgrid(x_int_diff,y_int_diff)
	
	#Mapping interpolation to interpolation grid
	phi_slice_int_diff=zeros((len(x_int_diff),len(y_int_diff)))
	
	print
	print "Interpolating over first slice restricted to image:"
	
	int_err=0
	int_tot=0
	img_tot=0
	
	for i in range(len(x_int_diff)):
		for j in range(len(y_int_diff)):
			
			phi_slice_int_diff[i,j]=phi_int(x_int_diff[i],y_int_diff[j])
			
			#Collecting pixels in embryo
			if sqrt((i-center[0])**2+(j-center[1])**2)<radius_data:
				int_err=int_err+(phi_slice_int_diff[i,j]-im_reg_ICs[i,j])
				img_tot=img_tot+im_reg_ICs[i,j]
				int_tot=int_tot+phi_slice_int_diff[i,j]
				
		#Print Progress
		if signal==None:
			curr_perc=int(100*i/float(shape(x_int_diff)[0]))
			sys.stdout.write("\r%d%%" %curr_perc)  
			sys.stdout.flush()
	print
	
	#Calculate difference between interpolation and picture
	int_diff=im_reg_ICs-phi_slice_int_diff	
	int_diff=absolute(int_diff)
	
	#Maximum value of img
	max_img=amax(im_reg_ICs)
	min_img=amin(im_reg_ICs)
	
	#Plotting
	fig=plt.figure()
	conlevels=linspace(phi_int_min,phi_int_max,50)
	
	#Data Image 
	ax_im=fig.add_subplot(211)
	plt_img=ax_im.contourf(im_reg_ICs,levels=conlevels,vmax=max_img,vmin=min_img)
	mark_radius_data=ptc.Circle((center[0],center[1]),radius_data,fill=False,linewidth=3,color='r')
	ax_im.add_patch(mark_radius_data)
	fig.colorbar(plt_img,ticks=[min_img,(max_img-min_img)/2,max_img], orientation='horizontal')
	
	#Regulated Interpolation
	ax_int=fig.add_subplot(212)
	plt_int=ax_int.contourf(x_int_grid,y_int_grid,phi_slice_int,levels=conlevels,vmax=phi_int_max,vmin=phi_int_min)
	fig.colorbar(plt_int, ticks=[phi_int_min,(phi_int_max-phi_int_min)/2,phi_int_max], orientation='horizontal')
	mark_radius_data2=ptc.Circle((center[0],center[1]),radius_data,fill=False,linewidth=3,color='r')
	ax_int.add_patch(mark_radius_data2)
	#ax_int.scatter(x_slice,y_slice,c='k',marker='o')
	
	
	#Contour plot for differences between im_reg_ICs and interpolation
	#ax_diff=fig.add_subplot(233)
	#plt_diff=ax_diff.contourf(x_int_diff_grid,y_int_diff_grid,int_diff,cmap=cm.winter)
	#fig.colorbar(plt_diff, orientation='horizontal')
	#mark_radius_data3=ptc.Circle((center[0],center[1]),radius_data,fill=False,linewidth=3,color='r')
	#ax_diff.add_patch(mark_radius_data3)
	
	
	
	
	#---------------------------------------------
	#IC check for slice -50
	#---------------------------------------------
	
	#slice_ind=slice_inds[1]
	
	##Getting slice enties of x,y,phi
	#phi_slice=phi.value[slice_ind]
	#x_slice=x[slice_ind]
	#y_slice=y[slice_ind]
	
	##Interpolation
	#phi_int=interp.Rbf(x_slice, y_slice, phi_slice, function='linear')
	
	##Creating interpolation grid
	#x_int = linspace(min(x_slice), max(x_slice), 600)
	#y_int = linspace(min(y_slice), max(y_slice), 600)
	#x_int_grid,y_int_grid=meshgrid(x_int,y_int)
	
	##Mapping interpolation to interpolation grid
	#phi_slice_int=zeros((len(x_int),len(y_int)))
	
	#phi_int_max=0.
	#phi_int_min=inf
	
	#print
	#print "Interpolating over first slice:"
	
	#for i in range(len(x_int)):
		#for j in range(len(y_int)):
			
			#phi_slice_int[i,j]=phi_int(x_int[i],y_int[j])
			
			#if phi_int_min>phi_slice_int[i,j]:
				#phi_int_min=phi_slice_int[i,j]
			#elif phi_int_max<phi_slice_int[i,j]:
				#phi_int_max=phi_slice_int[i,j]
		
				
		
		##Print Progress
		#curr_perc=int(100*i/float(shape(x_int)[0]))
		#sys.stdout.write("\r%d%%" %curr_perc)  
		#sys.stdout.flush()

	
	
	##Interpolation on same grid as img for difference calculation
	#x_int_diff=linspace(1, 512, 512)
	#y_int_diff=linspace(1, 512, 512)
	
	#x_int_diff_grid,y_int_diff_grid=meshgrid(x_int_diff,y_int_diff)
	
	##Mapping interpolation to interpolation grid
	#phi_slice_int_diff=zeros((len(x_int_diff),len(y_int_diff)))
	
	#print
	#print "Interpolating over first slice restricted to image:"
	
	#for i in range(len(x_int_diff)):
		#for j in range(len(y_int_diff)):
			
			#phi_slice_int_diff[i,j]=phi_int(x_int_diff[i],y_int_diff[j])
						
		##Print Progress
		#curr_perc=int(100*i/float(shape(x_int_diff)[0]))
		#sys.stdout.write("\r%d%%" %curr_perc)  
		#sys.stdout.flush()
	
	##Calculate difference between interpolation and picture
	#int_diff=im_reg_ICs-phi_slice_int_diff
	
	#int_diff=absolute(int_diff)
	
	
	##Plotting
	
	#conlevels=linspace(phi_int_min,phi_int_max,50)
	
	##Data Image 
	#ax_im2=fig.add_subplot(234)
	#plt_img2=ax_im2.contourf(im_reg_ICs,levels=conlevels,vmax=max_img,vmin=min_img)
	
	#fig.colorbar(plt_img2,ticks=[min_img,(max_img-min_img)/2,max_img], orientation='horizontal')
	#mark_radius_data4=ptc.Circle((center[0],center[1]),radius_data,fill=False,linewidth=3,color='r')
	#ax_im2.add_patch(mark_radius_data4)

	##Regulated Interpolation
	#ax_int2=fig.add_subplot(235)
	#plt_int2=ax_int2.contourf(x_int_grid,y_int_grid,phi_slice_int,levels=conlevels,vmax=phi_int_max,vmin=phi_int_min)
	#fig.colorbar(plt_int2,ticks=[phi_int_min,(phi_int_max-phi_int_min)/2], orientation='horizontal')
	#mark_radius_data5=ptc.Circle((center[0],center[1]),radius_data,fill=False,linewidth=3,color='r')
	#ax_int2.add_patch(mark_radius_data5)
	##ax_int2.scatter(x_slice,y_slice,c='k',marker='o')

	##Contour plot for differences between im_reg_ICs and interpolation
	#ax_diff2=fig.add_subplot(236)
	#plt_diff2=ax_diff2.contourf(x_int_diff_grid,y_int_diff_grid,int_diff,cmap=cm.winter)
	#fig.colorbar(plt_diff2, orientation='horizontal')
	#mark_radius_data6=ptc.Circle((center[0],center[1]),radius_data,fill=False,linewidth=3,color='r')
	#ax_diff2.add_patch(mark_radius_data6)
	
	plt.show()
	
	#fig.savefig("~/Documents/Research/PyFRAP/Code/simresults/blub.pdf",format='pdf', dpi=450)	
	
	print "int_err = ", int_err
	print "int_tot = ", int_tot
	print "img_tot = ", img_tot
	print "int_err/img_tot = ", int_err/img_tot
	raw_input()	
	return	

	