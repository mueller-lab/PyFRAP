#Script to analyze zstack images and compute geomerty

#=========================================================================================================================
#Importing modules
#=========================================================================================================================

from pyfrp_img_module import *
from pyfrp_misc_module import *
from pyfrp_gmsh_module import *
import skimage.io as skiio
import skimage.measure as skimsr
import matplotlib.pyplot as plt
import cv2
from scipy import ndimage
import os
import numpy as np

#=========================================================================================================================
#Parameters
#=========================================================================================================================

#~~~~~~~~~~~~~~~~~~
#Files

#Get current zstack
folder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/20160119_zstack_nobeads/01"
files=getSortedFileList(folder,".tif")
#files=files[8:]
outfile="zStack2"
outpath="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/outfiles/"

#Some flags what to do
debug_opt=1
debug_fill=1
sum_plot=1
apply_cutoff=0
apply_otsu=1
reverse_order=1
compose_polys=1

#Some empty vectors to save results
mpolys=[]

#Some parameters
stack_distance=4.51
volSize_px=25
fill_samples=2000
fill_mode="random"
dpoly=1

#Reconstruction algorithm
recon_alg="poisson"
#recon_alg="vcg"


meshlabscript="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/scripts/normals_"+recon_alg+".mlx"

#=========================================================================================================================
#Automatically Detecting Geometry Boundaries
#=========================================================================================================================


if sum_plot==1:
	#Creating figure for final plot
	fig_sum=plt.figure()
	fig_sum.show()
num_rows=ceil(floor(len(files))/4.)
j=1

#If reverse_order is selected, reverse order of files
if reverse_order==1:
	files.reverse()

files=files[7:12]

for fn in files:
		
	if debug_opt==1:
		#Creating plot window
		fig=plt.figure()
		fig.show()
		
	#Reading img file
	img_file=folder+"/"+fn
	img = cv2.imread(img_file)
	img=img[:,:,1]
	
	#Apply gaussian filter
	img=ndimage.gaussian_filter(img, sigma=2)
	
	if debug_opt==1:
		#Plotting img file
		ax=fig.add_subplot(241)
		ax.contourf(img)
	
	#Applying otsu algorithm
	if apply_otsu==1:
		threshhold,thresh=otsuImageJ(img,1,0,0)
	elif apply_cutoff==1:
		threshhold=0.3*img.max()
		thresh=img
		
		thresh[(where(img<threshhold))]=0		
	else:
		threshhold=0.76*img.max()
		thresh=img
		
	if debug_opt==1:
		#Plotting otsu result file
		ax2=fig.add_subplot(242)
		ax2.contourf(thresh,cmap="Greys")
	
	print "---------------------------------"
	print "Img:", fn
	
	#Finding contours
	if debug_opt==1:
		#Plot img
		ax3=fig.add_subplot(243)
		ax3.contourf(img,cmap="Greys")
		ax4=fig.add_subplot(244)
		ax4.contourf(img,cmap="Greys")
	
	kernel = np.ones((10, 10), np.uint8)
	thresh = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel)
	
	thresh = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE,kernel, iterations=4)
	contours=skimsr.find_contours(thresh,0)
	
	#Interpolate contours
	for i,cont in enumerate(contours):
		
	
	
	
	
	polys=[]
	for i,cont in enumerate(contours):
		poly = skimsr.approximate_polygon(cont, tolerance=dpoly)
		polys.append(poly)
		
		
		#Plot contours
		if debug_opt==1:
			ax3.plot(cont[:, 1], cont[:, 0],'g-', linewidth=2)
			ax4.plot(poly[:, 1], poly[:, 0],'r-', linewidth=2)
			
			ax4.plot(cont[0,1], cont[0,0],'b*', linewidth=2)
			ax4.plot(cont[-1,1], cont[-1,0],'m*', linewidth=2)
			
	
	print len(polys)
	
	#find longest polygon
	if not compose_polys:
		lens=[]
		for poly in polys:
			lens.append(len(poly[:,1]))
		
		ind=lens.index(max(lens))
		mpoly=polys[ind]
	
	else:
		if len(polys)>1:
	
			polysLeft=list(polys[1:])
	
			#End polygon
			mpoly=[]
			
			#Take first polygon to start with
			poly=polys[0]
			
			#Choose most south point
			pStart=poly[np.where(poly[:,1]==min(poly[:,1])),:]
			pStart=pStart[0][0]
			
			#Compute distance from pStart to all other points in contour
			d=np.linalg.norm(np.asarray(pStart)-poly,axis=1)
	
			#Sort according to distance to pStart
			sd=np.argsort(d)
			poly=poly[sd]

			#Define last point
			pEnd=poly[-1,:]
			
			#Define as mpoly
			mpoly=poly.copy()
			
			#Debugging plots
			if debug_opt==1:
				#ax4.plot(pStart[0], pStart[1],'b*', linewidth=2)
				#ax4.plot(pEnd[0], pEnd[1],'m*', linewidth=2)
				pass
			
			#Now loop through polysLeft and find the one that is closest
			while len(polysLeft)>0:
				
				distPStart=[]
				for poly in polysLeft:
					
					#Choose the point that is closest to the last point of previous contour
					d=np.linalg.norm(np.asarray(pEnd)-poly,axis=1) 
					
					#Sort according to distance to pEnd
					sd=np.argsort(d)
					poly=poly[sd]
					
					#Choose closest point
					pStart=poly[0,:]
					
					#Calculate its distance and append
					#print pStart
					#print pEnd
					
					distPStart.append(np.linalg.norm(pEnd-pStart))
				
				#Take the poly with the closest startpoint
				idx=distPStart.index(min(distPStart))
				poly=polysLeft[idx]
				
				#Define first point of poly
				pStart=poly[0,:]
				
				#Compute distance from pStart to all other points in contour
				d=np.linalg.norm(np.asarray(pStart)-poly,axis=1)
		
				#Sort according to distance to pStart
				sd=np.argsort(d)
				poly=poly[sd]

				#Define last point
				pEnd=poly[-1,:]
				
				#Throw out poly of polysLeft
				polysLeft.pop(idx)
				
				#Concatenate
				mpoly=np.concatenate((mpoly,poly),axis=0)
						
		else:
			mpoly=poly
		
	mpolys.append(mpoly)
	
	if debug_opt==1:
		ax5=fig.add_subplot(245)
		ax5.contourf(img,cmap="Greys")
		ax5.plot(mpoly[:, 1], mpoly[:, 0],'r-', linewidth=2)	
		
		plt.draw()
		raw_input()
	
	if sum_plot==1:
		ax_sum=fig_sum.add_subplot(4,num_rows,j)
		ax_sum.contourf(img,cmap="Greys")
		ax_sum.plot(mpoly[:, 1], mpoly[:, 0],'r-', linewidth=2)	
		ax_sum.set_title(fn)
		
	j=j+1

if sum_plot==1:	
	plt.draw()
	raw_input("Done finding contours, press ENTER to continue")

#=========================================================================================================================
#Generating random points to close 1st and last stack
#=========================================================================================================================

inpoints_all=[]

if debug_fill==1:
	fig=plt.figure()
	fig.show()
	j=1

for i in range(len(files)):
	
	if i==0 or i==len(files)-1:
		inpoints=[]
		if fill_mode=="random":
		
			npts=fill_samples
			
			while len(inpoints)<npts:
			
				rx=random.random() 	
				ry=random.random() 	
				
				x=mpolys[i][:,1].min()+rx*(mpolys[i][:,1].max()-mpolys[i][:,1].min())
				y=mpolys[i][:,0].min()+ry*(mpolys[i][:,0].max()-mpolys[i][:,0].min())
				
				poly=zip(list(mpolys[i][:,1]),list(mpolys[i][:,0]))
				
				if point_inside_polygon(x,y,poly):
					
					inpoints.append([x,y])
		
		elif fill_mode=="regular":
			
			dx=abs(mpolys[i][:,1].min()-mpolys[i][:,1].max())
			dy=abs(mpolys[i][:,0].min()-mpolys[i][:,0].max())
			
			d=sqrt((dx*dy)/fill_samples)
				
			xvec=arange(mpolys[i][:,1].min(),mpolys[i][:,1].max(),d)
			yvec=arange(mpolys[i][:,0].min(),mpolys[i][:,0].max(),d)
			
			poly=zip(list(mpolys[i][:,1]),list(mpolys[i][:,0]))
				
			for x in xvec:
				for y in yvec:
					if point_inside_polygon(x,y,poly):
						
						inpoints.append([x,y])
			
		inpoints_all.append(inpoints)
		
		if debug_fill==1:
		
			ax=fig.add_subplot(2,1,j)
			ax.plot(mpolys[i][:, 1], mpolys[i][:, 0],'r-', linewidth=2)
			
			for pt in inpoints:
				ax.plot(pt[0],pt[1],'g*')
			j=j+1
	
if debug_fill==1:
	plt.draw()
	raw_input("Done filling first and last zstack, press ENTER to continue")
		
#=========================================================================================================================
#Writing point cloud into .ply file
#=========================================================================================================================


ns=[]
for i in range(shape(mpolys)[0]):
	ns.append(shape(mpolys[i])[0])
nvertex=int(sum(ns))+int(len(inpoints_all[0]))+int(len(inpoints_all[1]))


f = open(outpath+"ply/"+outfile+".ply", 'w')

f.write("""ply
format ascii 1.0
""")
f.write("element vertex " + str(nvertex) + "\n") 
f.write("""property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
end_header
""")
	
for i in range(len(files)):
	red=int(float(i)/float(len(files))*255)
	blue=int(float(i)/float(len(files))*255)
	green=int(255-float(i)/float(len(files))*255)
	
	if i==0 or i==len(files)-1:
		
		if i==0:
			inpoints=inpoints_all[0]
			
		else:
			inpoints=inpoints_all[1]
			
			
		for k in range(shape(inpoints)[0]):
				
			f.write(str(inpoints[k][0]).replace('.',',') + " " + str(inpoints[k][1]).replace('.',',') + " " + str(-i*stack_distance).replace('.',',')+ " " + str(red) + " " + str(green) + " " + str(blue)  )  
			f.write("\n")
	
	
	for k in range(shape(mpolys[i])[0]):
			
		f.write(str(mpolys[i][k,1]).replace('.',',') + " " + str(mpolys[i][k,0]).replace('.',',') + " " + str(-i*stack_distance).replace('.',',')+ " " + str(red) + " " + str(green) + " " + str(blue)  )  
		f.write("\n")
	

f.close()

#=========================================================================================================================
#Passing .ply file to meshlab and execute meshlab script
#=========================================================================================================================

meshlabcmd="meshlabserver -s " + meshlabscript + " -i " + outpath+"ply/"+outfile+".ply -o "  + outpath +"stl/"+outfile + ".stl"

os.system(meshlabcmd)

#=========================================================================================================================
#Merging .stl file with .geo script
#=========================================================================================================================

f = open(outpath+"geo/"+outfile+".geo", 'w')

f.write('Merge "' + outpath+"stl/"+outfile + '.stl"; \n')
f.write("""Surface Loop(2) = {1};
Volume(3) = {2};
"""
)

f.close()

#=========================================================================================================================
#Calling gmsh
#=========================================================================================================================

os.system("gmsh -3 " + outpath+"geo/"+outfile +".geo -o " + outpath+"msh/"+outfile + ".msh")

#os.system("gmsh -3 " + outpath + ".msh")


	