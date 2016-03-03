#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Obsolete functions module for image analysis and simulation for PyFRAP toolbox, including following functions:
#1) 
 
#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

def gen_reg_solution(embryo,phi):
	
	start_time_reg=time.clock()
	
	slice_rad=embryo.radius_embr_px
	
	#Only need to regularize if we plot
	if embryo.plot_cont==1 or embryo.plot_all==1:
	
		#Creating regularization grid
		x_grid=linspace(embryo.center_embr_px[0]-slice_rad,embryo.center_embr_px[0]+slice_rad,embryo.res_reg)
		y_grid=linspace(embryo.center_embr_px[1]-slice_rad,embryo.center_embr_px[1]+slice_rad,embryo.res_reg)
		
		phi_reg=[]
			
		#-------------------------------------------------------------------------------------------------------------------
		#Finding regularization mapping
		#-------------------------------------------------------------------------------------------------------------------
			
		embryo.mesh_maps=[]
		start_time_map=time.clock()
		for i in range(shape(embryo.slice_height_px)[0]):
				
			#Only find mapping if not already given
			embryo.mesh_maps.append(map_mesh(x_grid,y_grid,x_cut[i],y_cut[i]))
			
			print "Mapping of", shape(phi_cut[i])[0],  " points found in ", time.clock()-start_time_map
			print "Mesh mapping found after", time.clock()-start_time_total	
		
		#-------------------------------------------------------------------------------------------------------------------
		#Creating regularized solution using reg-mapping
		#-------------------------------------------------------------------------------------------------------------------

		x_grid,y_grid=meshgrid(x_grid,y_grid)
		for i in range(shape(embryo.slice_height_px)[0]):
			
			phi_reg.append(reg_sol(x_grid,y_grid,slice_rad,embryo.center_embr_px,phi_cut[i],embryo.mesh_maps[i],1,embryo.debug_reg))	
	
	print "Regularized", shape(phi_cut[i])[0],  " points in ", time.clock()-start_time_reg
	
	return embryo,phi_reg,x_grid,y_grid

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds mapping from irregular gmsh mesh into regular nparray meshgrid for plotting
#x_cut,y_cut: lists of irregular mesh coordinates
#x_grid, y_grid: numpy arrays for regular mesh
#returns mapping matrix with rows [entry in irregular list, x coordinate index, y coordinate index]

def map_mesh(x_grid,y_grid,x_cut,y_cut):
	
	#Creating empty maping
	mesh_map=zeros(((shape(x_cut)[0]),3))
	
	#Finding mappings for each grid cell
	for i in range(shape(x_grid)[0]-1):
		for j in range(shape(y_grid)[0]-1):
			for k in range(shape(x_cut)[0]):
				if x_grid[i]<=x_cut[k] and x_cut[k]<=x_grid[i+1]:
					if y_grid[j]<=y_cut[k] and y_cut[k]<=y_grid[j+1]:
						mesh_map[k,:]=[int(k),int(i),int(j)]
	
	return mesh_map

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Uses mapping matrix from map_mesh to map phi_cut into regular mesh
#radius: radius of cylinder
#phi_cut: list of irregular mesh solution
#x_grid, y_grid: numpy arrays for regular mesh
#returns regularized solution phi_reg

def reg_sol(x_grid,y_grid,radius,center,phi_cut,mesh_map,smooth_opt,debug_reg):
	
	#Mapping solution entries
	phi_reg_sum=zeros(shape(x_grid))
	phi_reg_num=zeros(shape(x_grid))
	
	for i in range(shape(phi_cut)[0]):
		
		phi_reg_sum[mesh_map[i,1],mesh_map[i,2]]=phi_reg_sum[mesh_map[i,1],mesh_map[i,2]]+phi_cut[i]
		phi_reg_num[mesh_map[i,1],mesh_map[i,2]]=phi_reg_num[mesh_map[i,1],mesh_map[i,2]]+1
	
	if debug_reg==1:
		fig=plt.figure()
		ax1=fig.add_subplot(121)
		ax2=fig.add_subplot(122)
		
		ax1.contourf(phi_reg_sum)
		ax2.contourf(phi_reg_num)
		
		plt.show()
		
	#Averaging over cells
	phi_reg=zeros(shape(x_grid))
	
	empt_entr=[]
	
	for i in range(shape(x_grid)[0]):
		for j in range(shape(y_grid)[0]):
			#Entering nan outside of circle
			if sqrt((x_grid[i,j]-center[0])**2+(y_grid[i,j]-center[0])**2)>radius:
				phi_reg[i,j]=float('nan')
				
			else:
				if phi_reg_num[i,j]==0:
					if debug_reg==1:
						print "ERROR: Empty cell in regularization at x=", x_grid[i,j], " and y= ", y_grid[i,j], "."
					
					
						if smooth_opt==1:	
							print "Choose either lower res or smaller volSize"
						else:
							print "Choose either lower res, smaller volSize or turn on smooth option"
					
					phi_reg[i,j]=float('nan')
					
					#If smooth option chosen, add new empty cell
					if smooth_opt==1:
						empt_entr.append([i,j])
						
				else:
					phi_reg[i,j]=float(phi_reg_sum[i,j])/float(phi_reg_num[i,j])
	smooth_opt=0
	if smooth_opt==1:
		if debug_reg==1:
			print "There were" ,shape(empt_entr)[0], "empty entries out of", shape(x_grid)[0]**2, "cells"
		phi_reg=smooth_reg_sol(empt_entr,phi_reg,debug_reg)
	
	return phi_reg

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Smoothes out empty cells in phi_reg
#empt_entr: matrix with all nan entries in phi_reg that are within radius 
#returns smoothed out regularized solution phi_reg

def smooth_reg_sol(empt_entr,phi_reg,debug_reg):
	smooth_range=1
	smooth_count=0
	for i in range(shape(empt_entr)[0]):
		
		res=shape(phi_reg)[0]

		empt_x=empt_entr[i][0]
		empt_y=empt_entr[i][1]
		
		#Adding up all nan-entries around current cell
		phi_surr_num=0
		phi_surr_sum=0
		
		nan_count=0
		
		for j in range(-smooth_range,smooth_range+1):
			for k in range(-smooth_range,smooth_range+1):
				#Check if surrounding cell is out of grid
				if empt_x+j<=res-1 and empt_x+j>=0 and empt_y+k<=res-1 and empt_y+k>=0: 
					if isnan(phi_reg[empt_x+j,empt_y+k]):
						#Do nothing b/c cell empty
						nan_count=nan_count+1
					else:
						phi_surr_sum=phi_surr_sum+phi_reg[empt_x+j,empt_y+k]
						phi_surr_num=phi_surr_num+1
			
		#Check if all surrounding cell were empty, if not take arithmetic average for smoothening
		if phi_surr_num==0:
			if debug_reg==1:
				print "ERROR: No smoothing out possible, try larger values for smooth_range!" 
		else:
			phi_reg[empt_x,empt_y]=float(phi_surr_sum)/float(phi_surr_num)
			smooth_count=smooth_count+1
	
	if debug_reg==1:		
		print "Smoothed out", smooth_count, "points out of ", shape(empt_entr)[0], "possible points." 	
		
	return phi_reg


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculates volume of a selected region by summing up over cell volumes

def calc_vol_mesh(mesh,region):
	
	#----------------------------------------------------------------------------------------------------------------------------------------------------------
	#mesh = fipy mesh object
	#region = region to be summed up. region can be either given as a list of indices of cell centers, by [[xmin,xmax],[ymin,ymax],[zmin,zmax]], or [] indicating whole mesh
	#----------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#Getting cellcenters
	x,y,z=mesh.cellCenters
	
	#========================================================
	#Collecting all indices in region
	#========================================================
	
	#Check if a region is given:
	if shape(region)[0]<1:
		print "You selected to compute the volume of the whole geometry returned in vol_out"
		in_reg_ind=[]
		out_reg_ind=range(shape(x)[0])
	else:
		#Check how region is given
		if isinstance(region[0],list):
			#region given as list of min/max coordinates
			
			#Need to find all vertices that lie in region
			#Entries of list can be ":" indicating no limit
			xmin=region[0][0]
			if xmin==":":
				xmin=min(x)	
			xmax=region[0][1]
			if xmax==":":
				xmax=max(x)	
			ymin=region[1][0]
			if ymin==":":
				ymin=min(y)	
			ymax=region[1][1]
			if ymax==":":
				ymax=max(y)
			zmin=region[2][0]
			if zmin==":":
				zmin=min(z)	
			zmax=region[2][1]
			if zmax==":":
				zmax=max(z)
			
			in_reg_ind=[]
			out_reg_ind=[]
			
			for i in range(shape(x)[0]):
				if xmin<=x[i] and x[i]<=xmax and ymin<=y[i] and y[i]<=ymax and zmin<=z[i] and z[i]<=zmax:
					in_reg_ind.append(i)
				else:
					out_reg_ind.append(i)
		else: 
			#highlighted region given as list of nodes
			in_reg_ind=region
			
			#Creating out_reg_ind by systemitically removing in_reg_ind from range
			out_reg_ind=range(shape(x)[0])
			for ind in in_reg_ind:
				out_reg_ind.remove(ind)	
	
	print "There are ", shape(in_reg_ind)[0], "nodes inside region and", shape(out_reg_ind)[0], "outside the region" 
	
	#========================================================
	#Finally calculating volume in and outside of region
	#========================================================

	
	vol_in=sum(mesh.cellVolumes[in_reg_ind])
	vol_out=sum(mesh.cellVolumes[out_reg_ind])
	
	return vol_in, vol_out

##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Getting indices inside and outside of square, slice and rim

def get_ind_regions(offset_sq,side_length_sq,radius,center,rim,add_rim_from_radius,res,debug=False):
	
	#Empty index vectors
	ind_circ_x=[]
	ind_circ_y=[]
	ind_sq_x=[]
	ind_sq_y=[]
	ind_rim_x=[]
	ind_rim_y=[]
	
	if debug:
		ind_rim_debug=zeros((res,res))
		ind_squ_debug=zeros((res,res))
		ind_out_debug=zeros((res,res))
		ind_slice_debug=zeros((res,res))
		
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in square
			if check_square(i+1,j+1,offset_sq,side_length_sq):
				ind_sq_x.append(i)
				ind_sq_y.append(j)
				
				if debug:
					ind_squ_debug[i,j]=1
					
			else:
				#Check if in circle
				if check_circle(i+1,j+1,center,radius):
					ind_circ_x.append(i)
					ind_circ_y.append(j)
					if debug:
						ind_out_debug[i,j]=1
				
				#Check if in rim
				if add_rim_from_radius==1:
					if check_circle(i,j,center,radius) and not check_circle(i,j,center,radius-rim):
						ind_rim_x.append(i)
						ind_rim_y.append(j)
						if debug:
							ind_rim_debug[i,j]=1
				
				elif add_rim_from_radius==0:
					
					if check_circle(i,j,center,radius) and not check_circle(i,j,center,res/2-rim):
						
						ind_rim_x.append(i)
						ind_rim_y.append(j)
						if debug:
							ind_rim_debug[i,j]=1
	
	#Stacking ind_circ and ind_sq and get ind_slice
	ind_slice_x=ind_circ_x+ind_sq_x
	ind_slice_y=ind_circ_y+ind_sq_y
	
	if debug:
		print "======= get_ind_regions debugging output ======="
		print "Square is centered:", check_squ_centered_from_ind(ind_sq_x,ind_sq_y,res)
		print "Square has right size:", check_squ_size(ind_sq_x,ind_sq_y,side_length_sq)
		
	if debug:
		ind_slice=ind_out_debug+ind_squ_debug
		
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([2,2],titles=["squ","out","slice","rim"],sup="get_ind_regions debugging output")
		
		#Plot
		axes[0].contourf(ind_squ_debug)
		axes[1].contourf(ind_out_debug)
		axes[2].contourf(ind_slice_debug)
		axes[3].contourf(ind_rim_debug)
		
		plt.draw()
			
	return ind_circ_x,ind_circ_y,ind_sq_x,ind_sq_y, ind_slice_x, ind_slice_y, ind_rim_x, ind_rim_y

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Find indices of mesh centers in slice

def get_slice_indices(embryo):
	
	#Grabbing cellCenters of mesh
	x,y,z=embryo.mesh.cellCenters
	
	#Image resolution
	#res_img=shape(embryo.im_reg_ICs)[0]
	res_img=embryo.data_res_px
	
	#Empty lists
	slice_inds=[]
	slice_inds_img=[]
	
	#Grab img indices
	img_ind=where((x<=embryo.center_embr_px[0]+res_img/2) & (x>=embryo.center_embr_px[0]-res_img/2) & (y<=embryo.center_embr_px[1]+res_img/2) & (y>=embryo.center_embr_px[1]-res_img/2))[0]
	
	#Loop through all selected slices
	for j in range(shape(embryo.slice_height_px)[0]):
		
		#Grab slice indices			
		if embryo.slice_bottom==1:
			slice_ind=where((z<embryo.slice_height_px[j]+2*embryo.slice_width_px) & (z > embryo.slice_height_px[j]))[0]
				
		elif embryo.slice_bottom==0:
			slice_ind=where((z<embryo.slice_height_px[j]+embryo.slice_width_px) & (z > embryo.slice_height_px[j]-embryo.slice_width_px))[0]
		
		#Check which are the matching indices
		slice_ind_img=pyfrp_misc.match_vals(slice_ind,img_ind)			
			

		#Check if nodes in slice are too sparse
		if len(slice_ind)<50:
			print "WARNING: There are", len(slice_ind)," points in slice at z=", embryo.slice_height_px[j]
			print "There are two possible reasons for this: Either the slice is out of the geometry domain or your slice is not fine enough. Run debugging mode or refine your mesh."
		
		#Append results
		slice_inds.append(slice_ind)
		slice_inds_img.append(slice_ind_img)
		
	return slice_inds,slice_inds_img	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Getting indices inside and outside of square, slice and rim

def get_ind_circle(radius,center,res,debug=False):
	
	#Empty index vectors
	ind_circ_x=[]
	ind_circ_y=[]
	
	if debug:
		ind_slice_debug=zeros((res,res))
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in circle
			if check_circle(i+1,j+1,center,radius):
				ind_circ_x.append(i)
				ind_circ_y.append(j)
				if debug:
					ind_circ_debug[i,j]=1
	
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["slice"],sup="get_ind_circle debugging output")
		axes[0].contourf(ind_circ_debug)
		
	return ind_circ_x,ind_circ_y



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Getting concentrations in and outside of square and in slice

def analyzeDataset(analysis,signal=None,emb_count=None,debug=False):
	
	#Empty lists
	conc_sq=[]
	conc_circ=[]
	conc_slice=[]
	
	if debug or signal==None:
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "Starting analyzing dataset for concentration profile and regulated mesh"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Get indices for squ, slice and circ and extended region
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Fix center and offset of embryo for quad_red
	if embryo.quad_red:
		embryo2quad(embryo)
	
	#Indices of pixels in picture lying in out(circ), squ(sq) and slice
	ind_circ_x,ind_circ_y,ind_sq_x,ind_sq_y,ind_slice_x,ind_slice_y, ind_rim_x,ind_rim_y=get_ind_regions(embryo.offset_bleached_px,embryo.side_length_bleached_px,embryo.radius_embr_px,embryo.center_embr_px,embryo.rim,embryo.add_rim_from_radius,embryo.data_res_px,debug=embryo.debug_analysis)
	
	#Find pixels that are outside of original images
	num_ext=find_pixels_ext(embryo.radius_embr_px,embryo.center_embr_px,embryo.data_res_px,embryo.add_rim_img,embryo.quad_red,debug=embryo.debug_analysis)
	
	#If quad_red is selected, find indices of pixels in first quadrant
	if embryo.quad_red:
		
		inds_circ_x,inds_circ_y=regions2quad([[ind_circ_x,ind_circ_y]],embryo.data_res_px,debug=embryo.debug_analysis)[0]
		inds_sq_x,inds_sq_y=regions2quad([[ind_sq_x,ind_sq_y]],embryo.data_res_px,debug=embryo.debug_analysis)[0]
		inds_slice_x,inds_slice_y=regions2quad([[ind_slice_x,ind_slice_y]],embryo.data_res_px,debug=embryo.debug_analysis)[0]
		inds_rim_x,inds_rim_y=regions2quad([[ind_rim_x,ind_rim_y]],embryo.data_res_px,debug=embryo.debug_analysis)[0]
		
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Loop through images and compute concentrations
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	for i in range(shape(embryo.file_list)[0]):
		
		#Compose filename
		fn_img=str(embryo.fn_datafolder+'/'+embryo.file_list[i])
			
		#Reading in image
		im_vals=load_img(fn_img,embryo.data_enc)
		
		#Check if skimage reads in image as 2D array, if not grab channel of image with maximum range
		if len(shape(im_vals))>2:
			im_vals,ind_max=get_max_range_channel(im_vals,debug=embryo.debug_analysis)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Process image / Get image ready for readout
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		#show_process_poss(im_vals,embryo.fn_preimage,embryo.data_enc)
		
		im_vals = process_img(im_vals,embryo.fn_preimage,embryo.data_enc,quad_red=embryo.quad_red,flip_before_process=embryo.flip_before_process,norm_by_pre=embryo.norm_by_pre,gaussian=embryo.gaussian,gaussian_sigma=embryo.gaussian_sigma,data_offset=embryo.data_offset,debug=embryo.debug_analysis)
		
		#median_filter(im_vals,debug=True)
		#median_filter(im_vals,radius=5,debug=True)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Compute concentrations
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		#Get concentration in rim
		conc_rim=mean_conc(ind_rim_x,ind_rim_y,im_vals)
		
		#Get concentration in square, circle and slice
		conc_sq.append(mean_conc(ind_sq_x,ind_sq_y,im_vals,debug=embryo.debug_analysis))
		conc_circ.append(mean_ext_conc(ind_circ_x,ind_circ_y,im_vals,conc_rim,num_ext,embryo.add_rim_img))
		conc_slice.append(mean_ext_conc(ind_slice_x,ind_slice_y,im_vals,conc_rim,num_ext,embryo.add_rim_img))
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Dumping some values already into embryo object
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		if i==0:
			#ICs in extra variable
			if embryo.quad_red:
				im_reg_ICs=conv_skio_to_np(flip_quad(im_vals))
			else:
				im_reg_ICs=conv_skio_to_np(im_vals)
				
			#Rim concentration to embryo
			embryo.conc_rim=conc_rim
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Print Progress
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		curr_perc=int(100*i/float(shape(embryo.file_list)[0]))
		
		if signal==None:
			sys.stdout.write("\r%d%%" %curr_perc)  
			sys.stdout.flush()
		else:	
			if emb_count==None:
				signal.emit(curr_perc)
			else:
				signal.emit(curr_perc,emb_count)
	
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Final debugging plots
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if embryo.debug_analysis:
		
		#Create figure
		if not embryo.add_rim_img:
			fig,axes = pyfrp_plt.make_subplot([1,2],titles=["Analyzed timeseries", "ICs"],sup="analayze_dataset debugging output")
		elif embryo.add_rim_img:
			fig,axes = pyfrp_plt.make_subplot([1,3],titles=["Analyzed timeseries", "ICs",],sup="analayze_dataset debugging output")
		
			#Extended image
			plt_ext=axes[2].contourf(im_ext)
			fig.colorbar(plt_ext, orientation='horizontal')
		
		#Plot concentration timeseries
		axes[0].plot(embryo.tvec_data,conc_sq,'b',label='squ')
		axes[0].plot(embryo.tvec_data,conc_circ,'r',label='out')
		axes[0].plot(embryo.tvec_data,conc_slice,'g',label='slice')
		
		#Plot ICs
		plt_ICs=axes[1].contourf(im_reg_ICs)
		fig.colorbar(plt_ICs, orientation='horizontal')
		
		#Mark regions
		mark_radius_data=ptc.Circle((embryo.center_embr_px[0],embryo.center_embr_px[1]),embryo.radius_embr_px,fill=False,linewidth=3,color='r')
		axes[1].add_patch(mark_radius_data)
		mark_bleached_data=ptc.Rectangle((embryo.offset_bleached_px[0],embryo.offset_bleached_px[1]),embryo.side_length_bleached_px,embryo.side_length_bleached_px,fill=False,linewidth=3,color='r')
		axes[1].add_patch(mark_bleached_data)
		
		#Draw
		plt.draw()
		raw_input()
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Add results to embryo object and return it
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	embryo.squ_av_data_d=conc_sq
	embryo.out_av_data_d=conc_circ
	embryo.slice_av_data_d=conc_slice
	embryo.im_reg_ICs=im_reg_ICs
	
	return embryo


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds theoretical pixels that could be filled up with rim concentration

def find_pixels_ext(radius,center,res,add_rim_img,quad_red,debug=False):

	num_ext=0
	
	if add_rim_img:
	
		#Going through all pixels that are in extended image
		for i in range(int(ceil(center[0]-radius)),int(floor(center[0]+radius))):
			for j in range(int(ceil(center[1]-radius)),int(floor(center[1]+radius))):
				
				#Check if in slice
				if check_circle(i,j,center,radius):
					
					#Check if pixels are out of original picture
					if i<0 or i>res or j<0 or j>res:
						
						if quad_red:
							if check_quad(i,j):
								num_ext=num_ext+1
						else:
							num_ext=num_ext+1
							
	return num_ext

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fill up extended pixels with rim concentration

def genExtImg(ROIs,center,conc_rim,imgVals):
	
	#Get maximum range of extension
	idxX=[]
	idxY=[]
	for r in ROIs:
		idxX=idxX+r.imgIdxX
		idxY=idxY+r.imgIdxY
	Xmin=min(idxX)
	Ymin=min(idxY)
	Xmax=max(idxX)
	Ymax=max(idxY)
	
	
	
	#Getting resolution
	res=shape(imgVals)[0]
	
	#Resolution of extended image
	res_ext=floor(center[0]+radius)-ceil(center[0]-radius)
	
	#Offset between original and extended image
	offset_ext=int((res_ext-res)/2)
	
	#Going through all pixels that are in extended image
	for i in range(Xmin,Xmax+1):
		for j in range(Ymin,Ymax+1):
				
			#Check if pixels are out of original picture
			if i<0 or i>res or j<0 or j>res:
				im_ext[i+offset_ext,j+offset_ext]=conc_rim
			else:
				#Already enter original img into extended img
				im_ext[i+offset_ext,j+offset_ext]=img[i,j]	
	
	return im_ext

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Readjust embryo object such that it works with quad_red

def embryo2quad(embryo,debug=False):
	
	if debug:
		printNote("======= embryo2quad debugging output =======")
	
	if not check_circle_centered(embryo.geometry.center,embryo.dataResPx):
		center_old=embryo.geometry.center
		embryo.geometry.center=centerCircle(embryo.dataResPx)
		
		if debug:
			print "Adjusted center from: ", center_old , " to " , embryo.geometry.center
	
	if not check_squ_centered(embryo.offset_bleached_px,embryo.side_length_bleached_px,embryo.data_res_px):
		offset_old=embryo.offset_bleached_px
		embryo.offset_bleached_px=center_squ(embryo.side_length_bleached_px,embryo.data_res_px)
		
		if embryo.debug_analysis:
			print "Adjusted offset from: ", center_old , " to " , embryo.offset_bleached_px
		
	return embryo

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Apply Data Ics for Simulation from file read out

def apply_data_ics(embryo):
	
	res=int(embryo.data_res_px)

	if embryo.apply_data>2:
		apply_data=embryo.apply_data
		embryo.apply_data=1
	
	map_rad=-5.*ones([res,res])
	
	#ICs for regularized mesh in squares
	phi_IC_sq=embryo.im_reg_ICs
	
	phi_IC_rad=[]
	phi_IC_rad_mesh=[]
	phi_IC_bleached=[]
	
	#Number of embryo.radius_embr_px steps
	rad_steps=ceil(embryo.radius_embr_px/embryo.rad_step_px)
	rad_steps=int(rad_steps)
	
	if embryo.debug_preproc==1:
		print "embryo.rad_step_px=", embryo.rad_step_px
		print "rad_steps=", rad_steps
		
	#============================================
	#MODE=1: Approximation of ICs by radial levels
	#============================================

	if embryo.apply_data==1:
		
		#Some empty vectors
		phi_IC_rad=zeros(rad_steps,)
		phi_IC_rad_num=zeros(rad_steps,)
		
		out_idx_app=[]
		squ_idx_app=[]
		slice_idx_app=[]
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Creating circular regularized ICs
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		for i in range(res):	
			for j in range(res):
				
				#IDEA: Check starting from the outside if  r_k<[i,j]<=r_k-1
				#Going through all radial layers from the outside to the inside. If [i,j] is finally in layer, break, so we don't add [i,j] to more than one layer!!!
				for k in range(rad_steps):
					
					
					if embryo.debug_preproc==1:
					
						#For debugging purpose we keep track of the mappings
						if k==0:
							if sqrt((float(i)-embryo.center_embr_px[0])**2+(float(j)-embryo.center_embr_px[1])**2)<(rad_steps-k)*embryo.rad_step_px:
								
								#Check if outside of bleached area
								if embryo.offset_bleached_px[0]>i or i>embryo.offset_bleached_px[0]+embryo.side_length_bleached_px or embryo.offset_bleached_px[1]>j or j>embryo.offset_bleached_px[1]+embryo.side_length_bleached_px: 
									
									out_idx_app.append([i,j])
								
								else:
									squ_idx_app.append([i,j])
									
								
								slice_idx_app.append([i,j])
					
					
					if sqrt((float(i)-embryo.center_embr_px[0])**2+(float(j)-embryo.center_embr_px[1])**2)>=(rad_steps-(k+1))*embryo.rad_step_px and sqrt((float(i)-embryo.center_embr_px[0])**2+(float(j)-embryo.center_embr_px[1])**2)<(rad_steps-k)*embryo.rad_step_px:
							
						if embryo.debug_preproc==0:
							slice_idx_app.append([i,j])
							
						#Check if outside of bleached area
						if embryo.offset_bleached_px[0]>i or i>embryo.offset_bleached_px[0]+embryo.side_length_bleached_px or embryo.offset_bleached_px[1]>j or j>embryo.offset_bleached_px[1]+embryo.side_length_bleached_px: 
							
						
							#Instantly creating circular regularized ICs
							phi_IC_rad[k-1]=phi_IC_rad[k-1]+phi_IC_sq[i,j]
							phi_IC_rad_num[k-1]=phi_IC_rad_num[k-1]+1
							
							#Creating mapping for further operations
							map_rad[i,j]=k
								
							break			
						
						else:
							map_rad[i,j]=0
				
						
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Averaging over radial levels
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		for k in range(rad_steps):
			
			if phi_IC_rad_num[k]>0:
				
				phi_IC_rad[k]=float(phi_IC_rad[k])/float(phi_IC_rad_num[k])
				
				#phi_IC_rad[k]=float(phi_IC_rad[k])/float(phi_IC_rad_num[k])
	
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Make sure that radial levels are monotonically increasing at the boundary
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		smooth_bounds=0
		
		#Defining that boundary layers are 10% of all the layers
		if smooth_bounds==1:
			N_bound_layers=int(1./10.*rad_steps)
			
			if N_bound_layers<2 and rad_steps>2:
				N_bound_layers=2
				
			for k in range(N_bound_layers-1):
				
				if phi_IC_rad[k+1]>phi_IC_rad[k]:
					phi_IC_rad[k]=phi_IC_rad[k+1]
		
		
	#============================================
	#MODE=0: Ideal ICs: phi_out=const.
	#============================================

	elif embryo.apply_data==0:
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Going through all coordinates and check if within boundaries
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		phi_IC_rad=0.
		phi_IC_rad_num=0.
		
		out_idx_app=[]
		squ_idx_app=[]
		slice_idx_app=[]
		
		for i in range(res):
			for j in range(res):
				
				#Check if [i,j] is within embryo boundaries
				if sqrt((float(i)-embryo.center_embr_px[0])**2+(float(j)-embryo.center_embr_px[1])**2)<=embryo.radius_embr_px:
					
					#Check if outside of bleached area
					if embryo.offset_bleached_px[0]>i or i>embryo.offset_bleached_px[0]+embryo.side_length_bleached_px or embryo.offset_bleached_px[1]>j or j>embryo.offset_bleached_px[1]+embryo.side_length_bleached_px: 
						
						
						phi_IC_rad=phi_IC_rad+phi_IC_sq[i,j]
						phi_IC_rad_num=phi_IC_rad_num+1.
						map_rad[i,j]=1
						
						#Adding to out index
						out_idx_app.append([i,j])
			
					else:
						
						#Adding to square index
						squ_idx_app.append([i,j])
					
					#Adding to slice index
					slice_idx_app.append([i,j])
					
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Averaging 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		phi_IC_rad=phi_IC_rad/float(shape(slice_idx_app)[0])
	
	
	#============================================
	#MODE=2: Direct image application
	#============================================

	
	elif embryo.apply_data==2:
		
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Going through all coordinates and collect the ones that are in far outer ring
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		phi_IC_rad=0.
		phi_IC_rad_num=0.
		
		for i in range(res):
			for j in range(res):
				
				#Check if [i,j] is within embryo boundaries
				if sqrt((float(i)-embryo.center_embr_px[0])**2+(float(j)-embryo.center_embr_px[1])**2)<=embryo.radius_embr_px:
					
					#Check if [i,j] is outside 0.9*embryo.radius_embr_px to get outer layers
					if sqrt((float(i)-embryo.center_embr_px[0])**2+(float(j)-embryo.center_embr_px[1])**2)>0.9*embryo.radius_embr_px:
						
						
						phi_IC_rad=phi_IC_rad+phi_IC_sq[i,j]
						phi_IC_rad_num=phi_IC_rad_num+1.
						map_rad[i,j]=1
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Averaging 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		phi_IC_rad=phi_IC_rad/phi_IC_rad_num
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Finding entries in square
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	
	#Empty vectors
	map_bleached=zeros([res,res])
	phi_IC_bleached_num=0
	phi_IC_bleached_sum=0
	
	
	
	for i in range(res):	
		for j in range(res):
			if embryo.offset_bleached_px[0]<=i and i<=embryo.offset_bleached_px[0]+embryo.side_length_bleached_px and embryo.offset_bleached_px[1]<=j and j<=embryo.offset_bleached_px[1]+embryo.side_length_bleached_px: 
				map_bleached[i,j]=1
				phi_IC_bleached_sum=phi_IC_bleached_sum+phi_IC_sq[i,j]
				phi_IC_bleached_num=phi_IC_bleached_num+1
								
	phi_IC_bleached=float(phi_IC_bleached_sum)/float(phi_IC_bleached_num)
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Mapping back to mesh
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	
	phi_IC_rad_mesh=zeros([res,res])
	for i in range(res):
		for j in range(res):
			if map_bleached[i,j]==1:
				phi_IC_rad_mesh[i,j]=phi_IC_bleached
			
			elif map_rad[i,j]>0:
				if embryo.apply_data==1:
					phi_IC_rad_mesh[i,j]=phi_IC_rad[map_rad[i,j]]
				elif embryo.apply_data==0 or embryo.apply_data==2:
					phi_IC_rad_mesh[i,j]=phi_IC_rad
			
			
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Test plotting
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	plot_ic_opt=embryo.debug_preproc
	
	if plot_ic_opt==1:
		#Uncomment if plotting is necessary
		min_phi=min([min(phi_IC_sq[0]),min(phi_IC_sq[1])])
		min_phi=0
		max_phi=max([max(phi_IC_sq[0]),max(phi_IC_sq[1])])
		
		if embryo.apply_data==0:
			max_rad_phi=phi_IC_rad
			min_rad_phi=phi_IC_bleached
		else:
			max_rad_phi=max(phi_IC_rad)
			
			min_rad_phi=phi_IC_bleached
			
			print max_rad_phi, min_rad_phi
			raw_input()
			
		cflevels=linspace(min_phi,max_phi,rad_steps)
		cflevels_rad=linspace(min_rad_phi,max_rad_phi,rad_steps)
		
		
		x_grid=linspace(0,res,res)
		y_grid=linspace(0,res,res)
		x_grid,y_grid=meshgrid(x_grid,y_grid)
		
		fig_ICs=plt.figure()
		fig_ICs.show()
		ax_reg=fig_ICs.add_subplot(2,2,1)
		ax_rad=fig_ICs.add_subplot(2,2,2)
		ax_map_bl=fig_ICs.add_subplot(2,2,3)
		ax_map_rad=fig_ICs.add_subplot(2,2,4)
		
		
		reg_plot=ax_reg.contourf(x_grid,y_grid,phi_IC_sq,levels=cflevels)
		rad_plot=ax_rad.contourf(x_grid,y_grid,phi_IC_rad_mesh,levels=cflevels_rad)
		cbar = plt.colorbar(rad_plot)
		map_bl_plot=ax_map_bl.contourf(x_grid,y_grid,map_bleached)
		map_rad_plot=ax_map_rad.contourf(x_grid,y_grid,map_rad)
			
		plt.colorbar(map_rad_plot)
		plt.draw()
		raw_input()
	
	#Adding results to embryo object
	embryo.phi_IC_rad=phi_IC_rad
	embryo.phi_IC_rad_mesh=phi_IC_rad_mesh
	embryo.rad_steps=rad_steps
	embryo.phi_IC_bleached=phi_IC_bleached

	embryo.apply_data=apply_data
	
	return embryo	
	
	
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Use fipy to create frog mesh

def genFiPyBallMesh(volSize_px,radius_px,center):

	#Making some parameters local to pass it to fipy
	center_embr_px_x=center[0]
	center_embr_px_y=center[1]
	
	#Passing everything to fipy
	mesh = Gmsh3D(''' volSize_px = %(volSize_px)g;                                                                   r
		radius = %(radius_px)g;
		center_x = %(center_embr_px_x)g;
		center_y = %(center_embr_px_y)g;
		Point(1) = {center_x, center_y, -radius, volSize_px};
		Point(2) = {center_x+radius, center_y, -radius, volSize_px};
		Point(3) = {center_x, center_y+radius, -radius, volSize_px};
		Point(4) = {center_x-radius, cshrinkenter_y, -radius, volSize_px};
		Point(5) = {center_x, center_y-radius, -radius, volSize_px};
		Point(6) = {center_x, center_y, radius-radius, volSize_px};
		Point(7) = {center_x, center_y, -radius-radius, volSize_px};
		Circle(1) = {3, 1, 7};
		Circle(2) = {7, 1, 5};
		Circle(3) = {5, 1, 6};
		Circle(4) = {6, 1, 3};
		Circle(5) = {3, 1, 4};
		Circle(6) = {4, 1, 5};
		Circle(7) = {5, 1, 2};
		Circle(8) = {2, 1, 3};
		Circle(9) = {2, 1, 7};
		Circle(10) = {2, 1, 6};
		Circle(11) = {6, 1, 4};
		Circle(12) = {4, 1, 7};
		Line Loop(14) = {4, -8, 10};
		Ruled Surface(14) = {14};
		Line Loop(16) = {8, 1, -9};
		Ruled Surface(16) = {16};
		Line Loop(18) = {1, -12, -5};
		Ruled Surface(18) = {18};
		Line Loop(20) = {11, -5, -4};
		Ruled Surface(20) = {20};
		Line Loop(22) = {11, 6, 3};
		Ruled Surface(22) = {22};
		Line Loop(24) = {3, -10, -7};
		Ruled Surface(24) = {24};
		Line Loop(26) = {9, 2, 7};
		Ruled Surface(26) = {26};
		Line Loop(28) = {6, -2, -12};
		Ruled Surface(28) = {28};
		Surface Loop(30) = {20, 22, 28, 26, 16, 14, 24, 18};
		Volume(30) = {30};
		''' % locals()) 
	
	return mesh
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Use fipy to create fish mesh

def genFiPyDomeMesh(volSize_px,fish_inradius_px,fish_outradius_px,fish_dist_px,center):
	
	#Making some parameters local to pass it to fipy
	center_embr_px_x=center[0]
	center_embr_px_y=center[1]
	
	#Calculating some additional parameters
	z_interc=-((fish_outradius_px**2-fish_inradius_px**2)/(2*fish_dist_px) - fish_outradius_px + fish_dist_px/2)
	x_interc=sqrt(fish_inradius_px**2-((fish_outradius_px**2-fish_inradius_px**2)/fish_dist_px**2))
	x_interc=z_interc
	
	#Passing everything to fipy
	mesh = Gmsh3D(''' volSize_px = %(volSize_px)g; 
		fish_inradius_px = %(fish_inradius_px)g;
		fish_outradius_px = %(fish_outradius_px)g;
		fish_dist_px = %(fish_dist_px)g;
		z_interc = %(z_interc)g;
		x_interc = %(x_interc)g;
		center_x = %(center_embr_px_x)g;
		center_y = %(center_embr_px_y)g;
		Point(1) = {center_x, center_y, -fish_outradius_px-fish_dist_px, volSize_px};
		Point(2) = {center_x, center_y, -fish_outradius_px-fish_dist_px+fish_inradius_px, volSize_px};
		Point(7) = {center_x, center_y, -z_interc, volSize_px};
		Point(8) = {center_x+x_interc, center_y, -z_interc, volSize_px};
		Point(9) = {center_x, center_y+x_interc, -z_interc, volSize_px};
		Point(10) = {center_x-x_interc, center_y, -z_interc, volSize_px};
		Point(11) = {center_x, center_y-x_interc, -z_interc, volSize_px};
		Point(29) = {center_x, center_y, 0, volSize_px};
		Circle(12) = {8, 1, 2};
		Circle(13) = {9, 1, 2};
		Circle(14) = {10, 1, 2};
		Circle(15) = {11, 1, 2};
		Circle(24) = {8, 7, 9};
		Circle(25) = {9, 7, 10};
		Circle(26) = {10, 7, 11};
		Circle(27) = {11, 7, 8};
		Circle(30) = {8, 7, 29};
		Circle(31) = {9, 7, 29};
		Circle(32) = {10, 7, 29};
		Circle(33) = {11, 7, 29};
		Line Loop(35) = {33, -30, -27};
		Ruled Surface(35) = {35};
		Line Loop(37) = {32, -33, -26};
		Ruled Surface(37) = {37};
		Line Loop(39) = {31, -32, -25};
		Ruled Surface(39) = {39};
		Line Loop(41) = {24, 31, -30};
		Ruled Surface(41) = {41};
		Line Loop(43) = {12, -13, -24};
		Ruled Surface(43) = {43};
		Line Loop(45) = {13, -14, -25};
		Ruled Surface(45) = {45};
		Line Loop(47) = {14, -15, -26};
		Ruled Surface(47) = {47};
		Line Loop(49) = {15, -12, -27};
		Ruled Surface(49) = {49};
		Surface Loop(50) = {35, 37, 39, 41, 43, 49, 47, 45};
		Volume(51) = {50};
		''' % locals()) 
	
	return mesh
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Use fipy to create cylinder mesh

def genFiPyCylinderMesh(volSize_px,radius,height,center):
	
	#Making some parameters local to pass it to fipy
	center_embr_px_x=center[0]
	center_embr_px_y=center[1]
	
	#Passing everything to fipy
	mesh = Gmsh3D('''
		//Parameters
		volSize_px = %(volSize_px)g; 
		radius = %(radius)g;
		height = %(height)g;
		center_embr_px_x = %(center_embr_px_x)g;
		center_embr_px_y = %(center_embr_px_y)g;
		
		//Points of upper circle
		Point(1) = {center_embr_px_x, center_embr_px_y, 0, volSize_px}; 
		Point(2) = {-radius+center_embr_px_x, center_embr_px_y, 0, volSize_px};
		Point(3) = {center_embr_px_x, center_embr_px_y+radius, 0, volSize_px};
		Point(4) = {center_embr_px_x+radius, center_embr_px_y, 0, volSize_px};
		Point(5) = {center_embr_px_x, -radius+center_embr_px_y, 0, volSize_px};
		
		//Points of lower circle
		Point(6) = {center_embr_px_x, center_embr_px_y, -height, volSize_px}; 
		Point(7) = {-radius+center_embr_px_x, center_embr_px_y, -height, volSize_px};
		Point(8) = {center_embr_px_x, center_embr_px_y+radius, -height, volSize_px};
		Point(9) = {center_embr_px_x+radius, center_embr_px_y, -height, volSize_px};
		Point(10) = {center_embr_px_x, -radius+center_embr_px_y, -height, volSize_px};
		
		//Circles
		Circle(41) = {2, 1, 3}; 
		Circle(42) = {3, 1, 4}; 
		Circle(43) = {4, 1, 5}; 
		Circle(44) = {5, 1, 2};
		Circle(45) = {7, 6, 8}; 
		Circle(46) = {8, 6, 9}; 
		Circle(47) = {9, 6, 10}; 
		Circle(48) = {10, 6, 7};

		//Lines connecting circles 
		Line(49) = {2,7};
		Line(50) = {3,8};
		Line(51) = {4,9};
		Line(52) = {5,10};

		//Surfaces for disc
		Line Loop(53) = {44, 49, -48, -52};
		Ruled Surface(54) = {53};
		Line Loop(55) = {52, -47, -51, 43};
		Ruled Surface(56) = {55};
		Line Loop(57) = {42, 51, -46, -50};
		Ruled Surface(58) = {57};
		Line Loop(59) = {50, -45, -49, 41};
		Ruled Surface(60) = {59};
		Line Loop(61) = {42, 43, 44, 41};
		Ruled Surface(62) = {61};
		Line Loop(63) = {45, 46, 47, 48};
		Ruled Surface(64) = {63};
		Surface Loop(65) = {64, 60, 58, 62, 56, 54};

		//Volumes
		Volume(80) = {65};
		''' % locals()) 
	
	return mesh

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Use fipy to create fly mesh

def gen_fly_mesh(embryo):
	
	#Making some parameters local to pass it to fipy
	volSize_px=embryo.volSize_px
	ball_radius=embryo.fly_radius_px
	top_pt=embryo.fly_top_pt_px
	
	#NOTE: Need to insert formula for x_interc, z_interc here!
	
	center_embr_px_x=embryo.center_embr_px_x
	center_embr_px_y=embryo.center_embr_px_y
	
	embryo.mesh = Gmsh3D(''' volSize_px = %(volSize_px)g; 
		ball_radius = %(ball_radius)g;
		z_interc = %(z_interc)g;
		x_interc = %(x_interc)g;
		top_pt = %(top_pt)g;
		center_x = %(center_embr_px_x)g;
		center_y = %(center_embr_px_y)g;
		Point(1) = {center_x, center_y, 0, volSize_px};
		Point(2) = {center_x+ball_radius, center_y, -top_pt, volSize_px};
		Point(3) = {center_x, center_y+ball_radius, -top_pt, volSize_px};
		Point(4) = {center_x-ball_radius, center_y, -top_pt, volSize_px};
		Point(5) = {center_x, center_y-ball_radius, -top_pt, volSize_px};
		Point(6) = {center_x, center_y, -ball_radius-top_pt, volSize_px};
		Point(7) = {center_x, center_y, z_interc-top_pt, volSize_px};
		Point(8) = {cente Gmsh3r_x+x_interc, center_y, z_interc-top_pt, volSize_px};
		Point(9) = {center_x, center_y+x_interc, z_interc-top_pt, volSize_px};
		Point(10) = {center_x-x_interc, center_y, z_interc-top_pt, volSize_px};
		Point(11) = {center_x, center_y-x_interc, z_interc-top_pt, volSize_px};
		Point(29) = {center_x, center_y, top_pt-top_pt, volSize_px};
		Circle(12) = {2, 1, 3};
		Circle(13) = {3, 1, 4};
		Circle(14) = {4, 1, 5};
		Circle(15) = {5, 1, 2};
		Circle(16) = {2, 1, 6};
		Circle(17) = {3, 1, 6};
		Circle(18) = {4, 1, 6};
		Circle(19) = {5, 1, 6};
		Circle(20) = {2, 1, 8};
		Circle(21) = {3, 1, 9};
		Circle(22) = {4, 1, 10};
		Circle(23) = {5, 1, 11};
		Circle(24) = {8, 7, 9};
		Circle(25) = {9, 7, 10};
		Circle(26) = {10, 7, 11};
		Circle(27) = {11, 7, 8};
		Line(28) = {8, 29};
		Line(29) = {11, 29};
		Line(30) = {10, 29};
		Line(31) = {9, 29};
		Line Loop(33) = {19, -16, -15};
		Ruled Surface(33) = {33};
		Line Loop(35) = {18, -19, -14};
		Ruled Surface(35) = {35};
		Line Loop(37) = {18, -17, 13};
		Ruled Surface(37) = {37};
		Line Loop(39) = {17, -16, 12};
		Ruled Surface(39) = {39};
		Line Loop(41) = {14, 23, -26, -22};
		Ruled Surface(41) = {41};
		Line Loop(43) = {23, 27, -20, -15};
		Ruled Surface(43) = {43};
		Line Loop(45) = {20, 24, -21, -12};
		Ruled Surface(45) = {45};
		Line Loop(47) = {21, 25, -22, -13};
		Ruled Surface(47) = {47};
		Line Loop(49) = {26, 29, -30};
		Ruled Surface(49) = {49};
		Line Loop(51) = {29, -28, -27};
		Ruled Surface(51) = {51};
		Line Loop(53) = {24, 31, -28};
		Ruled Surface(53) = {53};
		Line Loop(55) = {31, -30, -25};
		Ruled Surface(55) = {55};
		Surface Loop(57) = {39, 37, 35, 33, 43, 41, 49, 51, 53, 45, 47, 55};
		Volume(57) = {57};
		''' % locals())

	return embryo



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if loops in .geo file are closed 

def check_ends(all_elmts,*curr_elmts):
	
	curr_elmts=list(curr_elmts)
	
	print "===================================================="
	
	#find curr_elmts in elmts
	elmts=[]
	i=1
	for el in curr_elmts:
		i=i+1 
		elmts.append(all_elmts[abs(el)])
		
	for i in range(1,len(elmts)):
		
		#Check if last entry of previous element is first entry of current element
		if (curr_elmts[i]<0 and elmts[i][0]!=elmts[i-1][0]) or (curr_elmts[i]>0 and elmts[i][0]!=elmts[i-1][-1]):
			
			print "Element ", elmts[i], "first entry does not match element", elmts[i-1], "last entry"
			
			#Revert signs
			if i==1:
				curr_elmts[i-1]=-curr_elmts[i-1]
				print "Element", curr_elmts[i-1], "has changed its sign"
				
			else:
				curr_elmts[i]=-curr_elmts[i]
				print "Element", curr_elmts[i], "has changed its sign"
		
	return curr_elmts

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out loop in .geo file

def printLoop(all_elmts,*curr_elmts):
	
	curr_elmts=list(curr_elmts)
	
	print "===================================================="
	
	#find curr_elmts in elmts
	elmts=[]
	for el in curr_elmts:
		if el<0:
			print "Element ID:", el, " is:", all_elmts[abs(el)][::-1]
		else:
			print "Element ID:", el, " is:", all_elmts[abs(el)]
	
	return

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prints out .geo-elements in and their ID

def printElements(elmts):
	
	for i in range(len(elmts)):
		
		print "Element ID:", i, "Element", elmts[i]
		
	return	