#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Plotting module for image analysis and simulation for PyFRAP toolbox, including following functions:
#(1)  geom_wire_plot: Produces wireframe plot of embryo geometry.
#(2)  surf_plot: Produces surface plot of given slice values.
#(3)  cont_plot: Produces contour plot of given slice values.
#(4)  conc_plot: Plots concentration profiles until step=step.
#(5)  all_plot: Produces figure with surf_plot, cont_plot and conc_plot as subfigures
#(6)  plot3Dmesh: Plot mesh in given region. Can also highlight given regions of the mesh. (slow)
#(7)  plot3Dtet: Plots a list of given tetrahedrons from a list of mesh nodes.
#(8)  plot3Dnodes: Plots nodes in given region similarly to plot3Dmesh.
#(9)  gen_masks: Generates masks of slice, square and out for plotting. Fills values out of region with fill_value (preferably 0 or NaN).

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#numpy
import numpy as np

#Plotting
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.collections as collections
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.tri as tri
#import pylab as plab

#Misc
import sys
import os

#PyFRAP
import pyfrp_img_module
from pyfrp_term_module import *

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Produces wireframe plot of embryo geometry

def geom_wire_plot(embryo,ax=None,gui=None):
	
	res_wire=40
	
	#~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~
	
	if gui!=None and ax!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	if ax==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111,projection='3d')
		
	elif gui!=None:
		ax=gui.ax
	
	#~~~~~~~~~~~~~~~~~
	#Make wireframe
	#~~~~~~~~~~~~~~~~~
	
	if embryo.geometry=="Cylinder":
		#Define Grid vectors
		x_wire=linspace(-embryo.cylinder_radius_px+embryo.center_embr_px[0]+0.0001,embryo.cylinder_radius_px+embryo.center_embr_px[0]-0.0001,res_wire)
		z_wire=linspace(0,-embryo.cylinder_height_px,res_wire)
		X_wire, Z_wire=meshgrid(x_wire,z_wire)
		Y_wire=embryo.center_embr_px[1]+sqrt(embryo.cylinder_radius_px**2-(X_wire-embryo.center_embr_px[0])**2)
		Y_wire2=-Y_wire+2*embryo.center_embr_px[1]
		
		ax.plot_wireframe(X_wire,Y_wire,Z_wire)
		ax.plot_wireframe(X_wire,Y_wire2,Z_wire)
		
	elif embryo.geometry=="Fish":
		
		#Define Grid vectors
		v_outer=math.acos(((embryo.fish_inradius_px**2-embryo.fish_outradius_px**2+embryo.fish_dist_px**2)/(2*embryo.fish_dist_px)-embryo.fish_dist_px)/embryo.fish_outradius_px)
		v_inner=math.acos(((embryo.fish_inradius_px**2-embryo.fish_outradius_px**2+embryo.fish_dist_px**2)/(2*embryo.fish_dist_px))/embryo.fish_inradius_px)
	
		#Parametric vectors inner
		u=r_[0:2*pi:50j]
		v=r_[v_inner:0:50j]
	
		x_wire_inner=embryo.center_embr_px[0]+embryo.fish_inradius_px*outer(cos(u),sin(v))
		y_wire_inner=embryo.center_embr_px[1]+embryo.fish_inradius_px*outer(sin(u),sin(v))
		z_wire_inner=-(embryo.fish_dist_px+embryo.fish_outradius_px)+embryo.fish_inradius_px*outer(ones(size(u)),cos(v))

		#Parametric vectors outer
		u=r_[0:2*pi:50j]
		v=r_[0:v_outer:50j]
		
		x_wire_outer=embryo.center_embr_px[0]+embryo.fish_outradius_px*outer(cos(u),sin(v))
		y_wire_outer=embryo.center_embr_px[1]+embryo.fish_outradius_px*outer(sin(u),sin(v))
		z_wire_outer=-embryo.fish_outradius_px+embryo.fish_outradius_px*outer(ones(size(u)),cos(v))
		
		ax.plot_wireframe(x_wire_outer,y_wire_outer,z_wire_outer,color='b')
		ax.plot_wireframe(x_wire_inner,y_wire_inner,z_wire_inner,color='b')
		
	elif embryo.geometry=="Frog":
		
		res_wire=20
		
		#Parametric vectors
		u=r_[0:2*pi:50j]
		v=r_[0:pi:50j]
		
		X_wire=embryo.center_embr_px[0]+embryo.frog_radius_px*outer(cos(u),sin(v))
		Y_wire=embryo.center_embr_px[1]+embryo.frog_radius_px*outer(sin(u),sin(v))
		Z_wire=-embryo.frog_radius_px+embryo.frog_radius_px*outer(ones(size(u)),cos(v))
		
		ax.plot_wireframe(X_wire,Y_wire,Z_wire,color='b')
	
	#Labels
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	
	#~~~~~~~~~~~~~~~~~
	#Draw
	#~~~~~~~~~~~~~~~~~
	
	if gui!=None:
		gui.canvas.draw()
	else:
		plt.draw()
		
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Produces surface plot of phi_cut at coordinates x_cut,y_cut

def surf_plot(embryo,x_cut,y_cut,phi_cut,min_phi=None,max_phi=None,ax=None,gui=None):
	
	#~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~
	
	if shape(embryo.slice_height_px)[0]>1:
		print "Warning: Surface plot is only available for a single slice."
		print "Will choose first slice for surface plot"
	
	#Compute max/min values of plot
	if min_phi==None:
		min_phi=min(phi_cut_surf)
		
	if max_phi==None:
		max_phi=max(phi_cut_surf)
	
	#Create figure and axis if necessary
	if gui!=None and ax!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	if ax==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111,projection='3d')
		#ax.set_zlim([min_phi,max()])
		
	elif gui!=None:
		ax=gui.ax	
	
	#Clearing axis
	ax.cla()
	
	ax.set_zlim([min_phi,max_phi])
	
	#~~~~~~~~~~~~~~~~~
	#Draw
	#~~~~~~~~~~~~~~~~~
	
	#surf_frame=ax.plot_trisurf(x_cut[0], y_cut[0], phi_cut[0], cmap=cm.jet, linewidth=0.2)
	surf_frame=ax.plot_trisurf(x_cut[0], y_cut[0], phi_cut[0], cmap=cm.jet, vmin=min_phi, vmax=max_phi, linewidth=0.2)
	
	if gui!=None:
		gui.canvas.draw()
	else:
		plt.draw()
		plt.pause(0.000000001)
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Produces contour plot of phi_cut at coordinates x_cut,y_cut

def cont_plot(embryo,x_cut,y_cut,phi_cut,cflevels=None,ax=None,gui=None):
	
	#~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~
	
	if shape(embryo.slice_height_px)[0]>1:
		print "Warning: Contour plot is only available for a single slice."
		print "Will choose first slice for surface plot"
	
	if gui!=None and ax!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	if ax==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
		
	elif gui!=None:
		ax=gui.ax	
	
	#Clearing axis
	ax.cla()
	
	#~~~~~~~~~~~~~~~~~
	#Draw
	#~~~~~~~~~~~~~~~~~
	
	cont_frame=ax.tricontourf(x_cut[0],y_cut[0],phi_cut[0],levels=cflevels)
	if gui!=None:
		gui.canvas.draw()
	else:
		plt.draw()
		plt.pause(0.000000001)
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Produces contour plot of phi_cut at coordinates x_cut,y_cut

def conc_plot(embryo,squ_av,out_av,slice_av,squ_small_av,out_small_av,step,ax=None,gui=None):
	
	#~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~
	
	if gui!=None and ax!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	if ax==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
		
	elif gui!=None:
		ax=gui.ax	

		
	#Creating time axis
	#tstep=range((step+2))
	#tstep = [tentry * timeStepDuration for tentry in tstep] 
	tstep=embryo.tvec_sim[:step+1]
	
	print shape(tstep), shape(squ_av)
	
	#~~~~~~~~~~~~~~~~~
	#Draw
	#~~~~~~~~~~~~~~~~~
	
	#Plotting
	squ_av_plt=ax.plot(tstep, squ_av,'b-',label="cin_large") 
	out_av_plt=ax.plot(tstep,out_av,'r-',label="cout_large")
	slice_av_plt=ax.plot(tstep,slice_av,'g-',label="cout_large")
	
	if embryo.avg_small==1:
		squ_av_small_plt=ax.plot(tstep,squ_small_av,'b--',label="cin_small")
		out_av_small_plt=ax.plot(tstep,out_small_av,'r--',label="cout_small")
		
	#Labels
	ax.set_ylabel('protein concentration')
	ax.set_xlabel('time')
	
	if gui!=None:
		gui.canvas.draw()
	else:
		plt.draw()
	
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Produces figure with surface, conc and cont plot
#NOTE: At some point might integrate wireframe in all_plot again, but for now this is not necessary

def all_plot(embryo,x_cut,y_cut,phi_cut,squ_av,out_av,slice_av,squ_small_av,out_small_av,step,axes=None,cflevels=None,max_phi=None,min_phi=None,ax=None,gui=None):
	
	#~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~
	
	if gui!=None and axes!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	
	if axes==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax_surf=fig.add_subplot(221,projection='3d')
		ax_cont=fig.add_subplot(222)
		ax_conc=fig.add_subplot(223)
	elif gui!=None:
		ax_surf=gui.ax_surf		
		ax_cont=gui.ax_cont	
		ax_conc=gui.ax_conc
	else:
		ax_surf=axes[0]
		ax_cont=axes[1]	
		ax_conc=axes[2]	
		
	if shape(embryo.slice_height_px)[0]>1:
		print "Warning: Surface plot is only available for a single slice."
		print "Will choose first slice for surface plot"
		
	ax_surf=surf_plot(embryo,x_cut,y_cut,phi_cut,min_phi=min_phi,max_phi=max_phi,ax=ax_surf,gui=None)
	ax_cont=cont_plot(embryo,x_cut,y_cut,phi_cut,cflevels=cflevels,ax=ax_cont,gui=None)
	ax_conc=conc_plot(embryo,squ_av,out_av,slice_av,squ_small_av,out_small_av,step,ax=ax_conc,gui=None)
	
	if gui!=None:
		gui.canvas.draw()
	else:
		plt.draw()
	
	return [ax_surf,ax_cont,ax_conc]

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creates 3D mesh plots for gmsh meshes

def plot3Dmesh(mesh,region,high_region,bound_region,alpha,markers,ax=None,gui=None):
	
	#----------------------------------------------------------------------------------------------------------------------------------------------------------
	#region options
	#1 -> plot whole mesh (blue)
	#2 -> plot only boundary vertices, can be restricted to bound_region, can be either given as a list of node indices or [[xmin,xmax],yxmin,ymax],[zmin,zmax]](green)
	#3 -> plot only highlighted region, high_region can be either given as a list of node indices or as a list [[xmin,xmax],yxmin,ymax],[zmin,zmax]] (red)
	#4 -> plot highlighted+whole mesh, mesh transparent
	#5 -> plot highlighted+boundary, boundary transparent
	#alpha = degree of transparency, 0=totally transparent, 1=not transparent
	#alpha = marker 
	#----------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if gui!=None and ax!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	if ax==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111, projection='3d')
		
	elif gui!=None:
		ax=gui.ax	
	
	#========================================================
	#Collecting all vertices to be plotted
	#========================================================
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Collecting all vertices in mesh
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if region==1:
		#Getting all Cell - to - Vertex IDs
		ordvert=mesh._getOrderedCellVertexIDs()
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	#Collecting all vertices in boundary that are in bound_region	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	elif region==2:
		#Getting all face vertices
		ordvert=mesh.faceVertexIDs
		#Check if a boundary region is given:
		if shape(bound_region)[0]<1:
			print "ERROR: No boundary region is given, going to plot whole boundary"
			in_bnd_ind=range(shape(ordvert)[1])
		else:
			#Check how region is given
			if isinstance(bound_region[0],list):
				#highlighted region given as list of min/max coordinates
				
				#Need to find all vertices that lie in highlighted region
				x,y,z=mesh.faceCenters
				
				#Entries of list can be ":" indicating no limit
				xmin=bound_region[0][0]
				if xmin==":":
					xmin=min(x)	
				xmax=bound_region[0][1]
				if xmax==":":
					xmax=max(x)	
				ymin=bound_region[1][0]
				if ymin==":":
					ymin=min(y)	
				ymax=bound_region[1][1]
				if ymax==":":
					ymax=max(y)
				zmin=bound_region[2][0]
				if zmin==":":
					zmin=min(z)	
				zmax=bound_region[2][1]
				if zmax==":":
					zmax=max(z)
				
				in_bnd_ind=[]
				out_bnd_ind=[]
				
				for i in range(shape(x)[0]):
					if xmin<=x[i] and x[i]<=xmax and ymin<=y[i] and y[i]<=ymax and zmin<=z[i] and z[i]<=zmax:
						in_bnd_ind.append(i)
					else:
						out_bnd_ind.append(i)
			else: 
				#highlighted region given as list of nodes
				in_bnd_ind=bound_region
				
				#Creating out_reg_ind by systemitically removing in_reg_ind from range
				out_bnd_ind=range(shape(ordvert)[1])
				for ind in in_bnd_ind:
					out_bnd_ind.remove(ind)	
		
		#Finally restricting vertices to bound_region
		ordvert=ordvert[:,in_bnd_ind]
	
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Collecting all vertices in high_region + all vertices
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	elif region==3 or region==4:
		#Getting all Cell - to - Vertex IDs
		ordvert=mesh._getOrderedCellVertexIDs()
		#Check if a highlighted region is given:
		if shape(high_region)[0]<1:
			print "ERROR: No highlighted region is given"
			raw_input()
		else:
			#Check how region is given
			if isinstance(high_region[0],list):
				#highlighted region given as list of min/max coordinates
				
				#Need to find all vertices that lie in highlighted region
				x,y,z=mesh.cellCenters
				
				#Entries of list can be ":" indicating no limit
				xmin=high_region[0][0]
				if xmin==":":
					xmin=min(x)	
				xmax=high_region[0][1]
				if xmax==":":
					xmax=max(x)	
				ymin=high_region[1][0]
				if ymin==":":
					ymin=min(y)	
				ymax=high_region[1][1]
				if ymax==":":
					ymax=max(y)
				zmin=high_region[2][0]
				if zmin==":":
					zmin=min(z)	
				zmax=high_region[2][1]
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
				in_reg_ind=high_region
				
				#Creating out_reg_ind by systemitically removing in_reg_ind from range
				out_reg_ind=range(shape(ordvert)[1])
				for ind in in_reg_ind:
					out_reg_ind.remove(ind)	
					
				
		
		print "There are", shape(in_reg_ind)[0] , " nodes in highlighted region and ", shape(out_reg_ind), "nodes outside of the region"
		
		if region==3:
			ordvert=ordvert[:,in_reg_ind]
		
		
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Collecting all vertices in high_region and bound_region
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	elif region==5:	
		ordvert=mesh.faceVertexIDs
		#Check if a boundary region is given:
		if shape(bound_region)[0]<1:
			print "ERROR: No boundary region is given, going to plot whole boundary"
			in_bnd_ind=range(shape(ordvert)[1])
		else:
			#Check how region is given
			if isinstance(bound_region[0],list):
				#highlighted region given as list of min/max coordinates
				
				#Need to find all vertices that lie in highlighted region
				x,y,z=mesh.faceCenters
				
				#Entries of list can be ":" indicating no limit
				xmin=bound_region[0][0]
				if xmin==":":
					xmin=min(x)	
				xmax=bound_region[0][1]
				if xmax==":":
					xmax=max(x)	
				ymin=bound_region[1][0]
				if ymin==":":
					ymin=min(y)	
				ymax=bound_region[1][1]
				if ymax==":":
					ymax=max(y)
				zmin=bound_region[2][0]
				if zmin==":":
					zmin=min(z)	
				zmax=bound_region[2][1]
				if zmax==":":
					zmax=max(z)
				
				in_bnd_ind=[]
				out_bnd_ind=[]
				
				for i in range(shape(x)[0]):
					if xmin<=x[i] and x[i]<=xmax and ymin<=y[i] and y[i]<=ymax and zmin<=z[i] and z[i]<=zmax:
						in_bnd_ind.append(i)
					else:
						out_bnd_ind.append(i)
			else: 
				#highlighted region given as list of nodes
				in_bnd_ind=bound_region
				
				#Creating out_reg_ind by systemitically removing in_reg_ind from range
				out_bnd_ind=range(shape(ordvert)[1])
				for ind in in_bnd_ind:
					out_bnd_ind.remove(ind)	
		
		ordvert=mesh._getOrderedCellVertexIDs()
		#Check if a highlighted region is given:
		if shape(high_region)[0]<1:
			print "ERROR: No highlighted region is given"
			raw_input()
		else:
			#Check how region is given
			if isinstance(high_region[0],list):
				#highlighted region given as list of min/max coordinates
				
				#Need to find all vertices that lie in highlighted region
				x,y,z=mesh.cellCenters
				
				#Entries of list can be ":" indicating no limit
				xmin=high_region[0][0]
				if xmin==":":
					xmin=min(x)	
				xmax=high_region[0][1]
				if xmax==":":
					xmax=max(x)	
				ymin=high_region[1][0]
				if ymin==":":
					ymin=min(y)	
				ymax=high_region[1][1]
				if ymax==":":
					ymax=max(y)
				zmin=high_region[2][0]
				if zmin==":":
					zmin=min(z)	
				zmax=high_region[2][1]
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
				in_reg_ind=high_region
				
				#Creating out_reg_ind by systemitically removing in_reg_ind from range
				out_reg_ind=range(shape(ordvert)[1])
				for ind in in_reg_ind:
					out_reg_ind.remove(ind)	
		
		#Getting ordered vertex IDs of highlighted regions
		ordvert=ordvert[:,in_reg_ind]
		
		#Getting vertex IDs from boundary
		facevert=mesh.faceVertexIDs
		
		#Creating fake IDs of -1 and attach to facevert so we can concatenate ordvert and facevert
		tempvec=-1*ones([1,shape(facevert)[1]])
		concatenate((facevert,tempvec),axis=0)
		facevert=concatenate((facevert,tempvec),axis=0)
		
		#Restricting facevert to whatever bound_region is:
		facevert=facevert[:,in_bnd_ind]
		
		ordvert=concatenate((ordvert,facevert),axis=1)
	
	#========================================================
	#Collecting all edges to be plotted
	#========================================================
	
	print "-------------------------------------------"
	print "Going through", shape(ordvert)[1], " Tetraeders"
	print "-------------------------------------------"
	
	edge_in_ind=[]
	edge_out_ind=[]
	
	#Extracting all IDs of vertices of edges to be plotted	
	edges=[]
	double_edges=0

	in_edges_removed=0
	out_edges_removed=0
	
	#Going through all tetraeders
	for j in range(shape(ordvert)[1]):
		curr_verts=list(ordvert[:,j])
		
		#Adding first vertex for closing of polygon
		curr_verts.append(curr_verts[0])
		
		#Going through all vetices of current tetraeders
		for i in range(shape(curr_verts)[0]-1):
			
			#Check if vertex ID is a fake ID used for concatenating
			if region==5:
				if curr_verts[i+1]==-1:
					break
			
			newedge=[curr_verts[i],curr_verts[i+1]]
			
			if newedge in edges:
				if region==3 or region==4:
					if j in in_reg_ind:
						in_edges_removed=in_edges_removed+1
					else:
						out_edges_removed=out_edges_removed+1
						
				#Edge already exists, count it and do nothing
				double_edges=double_edges+1
					
			else:
				#Actual new edge, add to collection of edges
				edges.append(newedge)
				
				#If highlighted regions, keep track which edge belongs to out and in
				if region==3: 
					edge_in_ind.append(shape(edges)[0]-1)
				elif region==4:
					if j in in_reg_ind:
						edge_in_ind.append(shape(edges)[0]-1)
					elif j in out_reg_ind:	
						edge_out_ind.append(shape(edges)[0]-1)
					else: 
						print "Something went wrone"
				elif region==5:
					#The first shape(in_reg_ind)[0] entries belong to highlighted region by the way ordvert was composed for region==5 
					if j<shape(in_reg_ind)[0]:
						edge_in_ind.append(shape(edges)[0]-1)
					else:
						#Rest belongs to the boundary vertices
						edge_out_ind.append(shape(edges)[0]-1)
						
		#Print Progress
		curr_perc=int(100*j/float(shape(ordvert)[1]))
		sys.stdout.write("\r%d%%" %curr_perc)  
		sys.stdout.flush()
		
	print "Removed", double_edges, " doubled edges from plotting"
	if region==3 or region==4:
		print "out of which", in_edges_removed, " were inside highlighted region and", out_edges_removed, "were outside regions" 
	
	#========================================================
	#Plotting all edges and vertices
	#========================================================
	
	print "-------------------------------------------"
	print "Started plotting"
	print "-------------------------------------------"
	print  shape(edges)[0], " to be plotted"
	
	for i in range(shape(edges)[0]):
		curr_edge=edges[i]
		
		#Creating vectors of coordinates of current edge
		oldvertx,oldverty,oldvertz=mesh.vertexCoords[:,curr_edge[0]]
		vertx,verty,vertz=mesh.vertexCoords[:,curr_edge[1]]
		xvec=[oldvertx,vertx]
		yvec=[oldverty,verty]
		zvec=[oldvertz,vertz]
		
		if region==1:
			style='b-'
			curr_alpha=alpha
		elif region==2:
			style='g-'
			curr_alpha=alpha
		elif region==3 or region==4:
			if i in edge_in_ind:
				style='r-'
				curr_alpha=1
			elif i in edge_out_ind:
				style='b-'
				curr_alpha=alpha
		elif region==5:		
			if i in edge_in_ind:
				style='r-'
				curr_alpha=1
			elif i in edge_out_ind:
				style='g-'
				curr_alpha=alpha
				
		ax.plot(xvec,yvec,zvec,style,alpha=curr_alpha,marker=markers)
		
		curr_perc=int(100*i/float(shape(edges)[0]))
		sys.stdout.write("\r%d%%" %curr_perc)  
		sys.stdout.flush()
	
	plt.show()
	
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creates node and mesh plot of a list of tetrahedrons

def plot3Dtet(mesh,high_region,ax=None,gui=None):
	
	#-----------------------------
	#NOTE: Tetrahedrons can only be given as list, not as geometrical region 
	
	ordvert=mesh._getOrderedCellVertexIDs()
	ordvert=ordvert[:,high_region]
	
	edges=[]
	double_edges=0
	
	#Going through all tetraeders
	for j in range(shape(ordvert)[1]):
		
		curr_verts=list(ordvert[:,j])
		
		
		#Going through all vetices of current tetraeders
		for i in range(shape(curr_verts)[0]-1):
			for k in range(i+1,shape(curr_verts)[0]):
				newedge=[curr_verts[i],curr_verts[k]]
				
				if newedge in edges:
							
					#Edge already exists, count it and do nothing
					double_edges=double_edges+1
						
				else:
					#Actual new edge, add to collection of edges
					edges.append(newedge)
		

	for i in range(shape(edges)[0]):
		curr_edge=edges[i]
		
		#Creating vectors of coordinates of current edge
		oldvertx,oldverty,oldvertz=mesh.vertexCoords[:,curr_edge[0]]
		vertx,verty,vertz=mesh.vertexCoords[:,curr_edge[1]]
		xvec=[oldvertx,vertx]
		yvec=[oldverty,verty]
		zvec=[oldvertz,vertz]
		
		if i==0:		
			ax.plot(xvec,yvec,zvec,'b-',label='Edge')
			ax.plot(xvec,yvec,zvec,'m*',label='Vertex')
			
		else:
			ax.plot(xvec,yvec,zvec,'b-')
			ax.plot(xvec,yvec,zvec,'m*')
			
		
	x,y,z=mesh.cellCenters()
	
	x=x[high_region]
	y=y[high_region]
	z=z[high_region]
	
	ax.plot(x,y,z,'ro',label='Node')
	ax.legend()
	plt.show()
				
	return ax


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creates 3D scatter plot of all nodes

def plot3Dnodes(mesh,region,high_region,ax=None,gui=None):
	
	#----------------------------
	#region options
	#1 -> plot all nodes (blue)
	#2 -> plot only boundary nodes (green)
	#3 -> plot only highlighted region, high_region can be either given as a list of node indices or as a list [[xmin,xmax],yxmin,ymax],[zmin,zmax]] (red)
	#4 -> plot highlighted+all nodes
	#5 -> plot highlighted + boundary
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Check inputs
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if gui!=None and ax!=None:
		print "Both ax and gui are not None, I am confused. Not going to plot"
		return
	if ax==None and gui==None:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111, projection='3d')
		
	elif gui!=None:
		ax=gui.ax	
	
	#========================================================
	#Collecting all nodes to be plotted
	#========================================================
	
	if region==1:
		#Getting all Cell - to - Vertex IDs
		x,y,z=mesh.cellCenters
	elif region==2:
		#Getting all face vertices
		x,y,z=mesh.getFaceCenters()	
	elif region==3 or region==4 or region==5:
		
		#Check if a highlighted region is given:
		if shape(high_region)[0]<1:
			print "ERROR: No highlighted region is given"
			raw_input()
		else:
			
			x,y,z=mesh.cellCenters
			#Check how region is given
			if isinstance(high_region[0],list):
				#highlighted region given as list of min/max coordinates
				
				#Need to find all vertices that lie in highlighted region
				
				
				#Entries of list can be ":" indicating no limit
				xmin=high_region[0][0]
				if xmin==":":
					xmin=min(x)	
				xmax=high_region[0][1]
				if xmax==":":
					xmax=max(x)	
				ymin=high_region[1][0]
				if ymin==":":
					ymin=min(y)	
				ymax=high_region[1][1]
				if ymax==":":
					ymax=max(y)
				zmin=high_region[2][0]
				if zmin==":":
					zmin=min(z)	
				zmax=high_region[2][1]
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
				in_reg_ind=high_region
				
				#Creating out_reg_ind by systemitically removing in_reg_ind from range
				out_reg_ind=range(shape(x)[0])
				for ind in in_reg_ind:
					out_reg_ind.remove(ind)	
			
		print "There are", shape(in_reg_ind)[0] , " nodes in highlighted region and ", shape(out_reg_ind), "nodes outside of the region"
		
		if region==3:
			
			x=x[in_reg_ind]
			y=y[in_reg_ind]
			z=z[in_reg_ind]
		
		elif region==4:
			xin=x[in_reg_ind]
			yin=y[in_reg_ind]
			zin=z[in_reg_ind]
	
			xout=x[out_reg_ind]
			yout=y[out_reg_ind]
			zout=z[out_reg_ind]
	
		elif region==5:
			xin=x[in_reg_ind]
			yin=y[in_reg_ind]
			zin=z[in_reg_ind]
			
			xout,yout,zout=mesh.getFaceCenters()
	
	#========================================================
	#Plotting nodes
	#========================================================
	
	print "-------------------------------------------"
	print "Started plotting"
	print "-------------------------------------------"
	print  shape(x)[0], "nodes to be plotted"
	
	if region==1:
		color='b'
	elif region==2: 
		color='g'
	elif region==3:
		color='r'
	
	if region<=3:
		ax.scatter(x, y, z,c=color)
	elif region==4:
		ax.scatter(xin, yin, zin,c='r',alpha=1)
		ax.scatter(xout, yout, zout,c='b',alpha=0.1)
	elif region==5:
		ax.scatter(xin, yin, zin,c='r',alpha=1)
		ax.scatter(xout, yout, zout,c='g',alpha=0.1)
	
	
	plt.show()
	
	return

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create mask for squ, slice, out, fills values outside of radius with fill_value

def gen_masks(grid_x,grid_y,fill_value,center,radius,offset,sidelength):

	mask_slice=zeros(shape(grid_x))
	mask_out=zeros(shape(grid_x))
	mask_squ=zeros(shape(grid_x))
	
	#Go through all interpolation coordinates
	for i in range(shape(grid_x)[0]-1):
		for j in range(shape(grid_y)[1]-1):
			#Check if in slice
			if sqrt((grid_x[i,j]-center[0])**2+(grid_y[i,j]-center[1])**2)<radius:
				mask_slice[i,j]=1
				#Check if in squ
				if offset[0]<=grid_x[i,j] and offset[1]<=grid_y[i,j] and grid_x[i,j]<=offset[0]+sidelength and grid_y[i,j]<=offset[1]+sidelength:
					mask_squ[i,j]=1	
				else:	
					mask_out[i,j]=1
			else:
				mask_out[i,j]=fill_value
				mask_squ[i,j]=fill_value
				mask_slice[i,j]=fill_value

	return mask_squ, mask_out, mask_slice

#===================================================================================================================================
#Opens figure with first post image and give click selector for boundaries
#===================================================================================================================================

class FRAPBoundarySelector():
	
	def __init__(self,embryo):
		
		#Passing embryo to class
		self.embryo=embryo	
		
		#Plot resolution
		self.dpi = 100
	
		#Some bookkeeping variables
		self.radiusPt=[]
		self.centerPt=[]
		self.centerMarker=[]
		
		#Results
		self.radius=None
		self.center=None
		
		#Creating figure and canvas
		self.createCanvas()
			
		#Plot image if existent
		self.showFirstDataImg()
		plt.show()
	
	def createCanvas(self):
			
		h=500/self.dpi
		v=500/self.dpi
		
		self.fig = plt.figure()
		self.fig.show()
		self.fig.set_size_inches(h,v,forward=True)
		
		self.canvas=self.fig.canvas
		
		self.ax = self.fig.add_subplot(111)
		self.ax.set_xlim([0, self.embryo.dataResPx])
		self.ax.set_ylim([0, self.embryo.dataResPx])
		
		#Connect to mouseclick event
		self.fig.canvas.mpl_connect('close_event', self.closeFigure)
		self.canvas.mpl_connect('button_press_event', self.getMouseInput)
		
		self.canvas.draw()
		
		return 
	
	def closeFigure(self,event):
		return self.center,self.radius
		
	def getEmbryo(self):	
		return self.embryo
	
	def getResults(self):
		return self.center,self.radius
	
	def drawCenterMarker(self):
		
		centerImg=[self.embryo.dataResPx/2.,self.embryo.dataResPx/2.]
		
		if len(self.centerMarker)>0:
			self.clearCenterMarker()
		else:	
			pt=ptc.Circle(centerImg,radius=3,fill=True,color='y')
			self.centerMarker.append(self.ax.add_patch(pt))
			
		self.fig.canvas.draw()
		
		return self.centerMarker
		
	def clearCenterMarker(self):
	
		for pt in self.centerMarker:
			pt.remove()
			
		self.fig.canvas.draw()	
			
		self.centerMarker=[]
	
	def drawCenter(self):
		
		if len(self.centerPt)>0:
			self.clearCenter()
	
		pt=ptc.Circle(self.center,radius=3,fill=True,color='r')
		self.centerPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.centerPt
	
	def clearCenter(self):
	
		for pt in self.centerPt:
			pt.remove()
			
		self.centerPt=[]
		
		self.fig.canvas.draw()
		
	def drawRadius(self):
		
		if len(self.radiusPt)>0:
			self.clearRadius()
			
		pt=ptc.Circle(self.center,radius=self.radius,fill=False,color='r')
		self.radiusPt.append(self.ax.add_patch(pt))
		
		self.fig.canvas.draw()
		
		return self.radiusPt
	
	def clearRadius(self):
	
		for pt in self.radiusPt:
			pt.remove()
			
		self.radiusPt=[]
		
		self.fig.canvas.draw()
		
		return self.radiusPt
		
		
	def showFirstDataImg(self):
		
		self.embryo.updateFileList()
		
		fnImg=self.embryo.fnDatafolder
		if fnImg[-1]!='/':
			fnImg=fnImg+'/'
		fnImg=fnImg+self.embryo.fileList[0]
		
		img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
	
		self.showImg(img)
	
	def showImg(self,img):
		
		self.ax.imshow(img)
		self.fig.canvas.draw()
		
		return self.canvas
	
	def computeRadiusFromCoordinate(self,x,y):
		return np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
		
	def getMouseInput(self,event):
		
		#Check if click in axes
		if event.xdata==None:
			return
		
		#Left click to define center
		if event.button==1:
			
			self.center=[event.xdata,event.ydata]
			
			self.drawCenter()
			
			if len(self.radiusPt)>0:
				self.drawRadius()
				
		#Right click to define radius
		elif event.button==3:
			
			if len(self.centerPt)>0:
			
				self.radius=self.computeRadiusFromCoordinate(event.xdata,event.ydata)
				self.drawRadius()
			
		#Middle click to activate center marker	
		if event.button==2:
			
			self.drawCenterMarker()
				
		self.fig.canvas.draw()
		
		return
	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Generates matplotlib figure with (x,y) subplots

def makeSubplot(size,titles=None,tight=False,sup=None,proj=None,fig=None):
	
	#How many axes need to be created
	n_ax=size[0]*size[1]
	
	#Check if titles has right size
	if titles!=None:
		
		if len(titles)!=n_ax:
			print "Warning: len(titles)!=n_ax, will not print titles!"
			titles=None
	if proj!=None:
		if len(proj)!=n_ax:
			print "Warning: len(proj)!=n_ax, will not do projections!"
			proj=n_ax*[None]
	else:
		proj=n_ax*[None]
	
	#Creating figure
	if fig==None:
		fig=plt.figure()
		fig.show()
	fig.set_tight_layout(tight)
	
	#Suptitle
	if sup!=None:
		fig.suptitle(sup)
	
	#Add axes
	axes=[]
	for i in range(n_ax):
		
		if proj[i]!=None:
			ax=fig.add_subplot(size[0],size[1],i+1,projection=proj[i])
		else:
			ax=fig.add_subplot(size[0],size[1],i+1)
		axes.append(ax)
		
		#Print titles 
		if titles!=None:
			ax.set_title(titles[i])
	
	#Draw
	plt.draw()
	
	#Return axes handle
	return fig,axes

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Adjust display range of imshow plots in list of axes

def adjustImshowRange(axes,vmin=None,vmax=None):
	
	#Loop through axes
	for ax in axes:
		
		#Grab plot
		implot=findImageArtist(ax,"Image")
		
		#Rescale
		if implot!=None:
			implot.set_clim(vmin,vmax)
		
	#Draw
	redraw(ax)
		
	return axes

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Finds artist with key in name and returns it

def findImageArtist(ax,key):
	
	c=ax.get_children()
	
	for x in c:
		if key in str(x):
			return x
			break
		
	return None	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Redraws axes's figure

def redraw(ax):
			
	ax.get_figure().canvas.draw()
	
	return ax

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Plot timeseries all-in-one function

def plotTS(xvec,yvec,label='',title='',sup='',ax=None,color=None,linewidth=1,legend=True,linestyle='-'):
		
	if len(xvec)!=len(yvec):
		printWarning('len(xvec) != len (yvec). This could be due to incomplete simulation/analysis/pinning. Will not plot.')
		return None
	
	if ax==None:
		fig,axes = makeSubplot([1,1],titles=[title],sup=sup,tight=False)
		ax=axes[0]
	
	ax.plot(xvec,yvec,color=color,label=label,linestyle=linestyle)
	
	if legend:
		ax.get_figure().subplots_adjust(right=0.7)
		ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	
	redraw(ax)
	
	return ax
		