#Script to analyze zstack images and compute geomerty

from pyfrp_img_module import *
from pyfrp_misc_module import *
from pyfrp_gmsh_module import *
import skimage.io as skiio
import matplotlib.pyplot as plt
import cv2
from scipy import ndimage
import os



#Get current zstack
folder="/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/Theresa/FRAP_10kDa-A488-Dextran_10Aga_20140520/exp43/geo"
files=get_sorted_folder_list(folder,".tif")

#files=files[:3]
#Some flags what to do
debug_opt=0
apply_thresh=1
reverse_order=1
force_larger=1

#Some empty vectors to save results
radiuses_final=[]
centers_final=[]

#Some parameters
stack_distance=5
volSize_px=25

#=========================================================================================================================
#Automatically Detecting Geometry Boundaries
#=========================================================================================================================


#Creating figure for final plot
fig_sum=plt.figure()
fig_sum.show()
num_rows=ceil(floor(len(files))/4.)
j=1

#If reverse_order is selected, reverse order of files
if reverse_order==1:
	files.reverse()

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
		ax=fig.add_subplot(131)
		ax.contourf(img)
	
	#Applying otsu algorithm
	if apply_thresh==1:
		threshhold,thresh=otsu_imagej(img,1,0,0)
	else:
		threshhold=0.76*img.max()
		thresh=img
		
	if debug_opt==1:
		#Plotting otsu result file
		ax2=fig.add_subplot(132)
		ax2.contourf(thresh,cmap="Greys")
	
	#Finding contours
	#contours=find_embr_bnds(img,img.max(),threshhold)
	contours,max_len,thresh_opt=find_longest_boundary(img,img.max())
	
	
	#com=compute_com(img)
	
	com=[256,256]
	
	radiuses=[]
	for cont in contours:
		for pt in cont:
			curr_r=sqrt((pt[0]-com[0])**2+(pt[1]-com[1])**2)
		
			if force_larger==1:
				if len(radiuses_final)>0:
					if curr_r>=0.9*radiuses_final[-1]:
						radiuses.append(curr_r)
				else:
					radiuses.append(curr_r)
			else:
				radiuses.append(curr_r)
	
	radius=compute_radius_hist(radiuses,0)
	
	radiuses_final.append(radius)
	centers_final.append(com)
	
	print "---------------------------------"
	print "Img:", fn
	print "Center:", com
	print "Radius:", radius
	
	
	if debug_opt==1:
		#Plot img
		ax3=fig.add_subplot(133)
		ax3.contourf(img,cmap="Greys")
		
		#Plot contours
		for cont in contours:
			ax3.plot(cont[:, 1], cont[:, 0],'g-', linewidth=2)
		
		#Plot final circle
		pt=ptc.Circle(com,radius=3,fill=True,color='r')
		ax3.add_patch(pt)
		circ=ptc.Circle(com,radius=radius,fill=False,color='r',linewidth=3)
		ax3.add_patch(circ)
		
		
		plt.draw()
		raw_input()
	
	
	ax_sum=fig_sum.add_subplot(4,num_rows,j)
	ax_sum.contourf(img,cmap="Greys")
	pt=ptc.Circle(com,radius=3,fill=True,color='r')
	ax_sum.add_patch(pt)
	circ=ptc.Circle(com,radius=radius,fill=False,color='r',linewidth=3)
	ax_sum.add_patch(circ)
	
	j=j+1
	
plt.draw()
raw_input("Done finding contours, press ENTER to continue")

#=========================================================================================================================
#Write everything into geo file
#=========================================================================================================================

f = open(folder+"/zStack.geo", 'w')

#Write general properties

f.write("volSize_px="+str(volSize_px)+";")
f.write("\n")

f.write("//Defining parameters of each stack")
f.write("\n")
#Write in height of stacks
for i in range(len(files)):
	
	f.write("//Stack number " + str(i))
	f.write("\n")
	f.write("radius_"+str(i)+"="+str(radiuses_final[i])+";")
	f.write("\n")
	f.write("center_x_"+str(i)+"="+str(centers_final[i][0])+";")
	f.write("\n")
	f.write("center_y_"+str(i)+"="+str(centers_final[i][1])+";")
	f.write("\n")
	f.write("height_"+str(i)+"="+str(i*stack_distance)+";")
	f.write("\n")

j=1
center_pts=[]
right_pts=[]
upper_pts=[]
left_pts=[]
lower_pts=[]
circles=[]

f.write("\n")
f.write("\n")

elements=[]
elements.append("Dummy")

#------------------
#Create circles
for i in range(len(files)):
	
	f.write("//Circle of stack " + str(i))
	f.write("\n")
	
	#~~~~~~~~~~~
	#Points
	#Center point
	f.write("Point("+str(j)+") = {center_x_"+str(i)+", center_y_"+str(i)+", height_"+str(i)+", volSize_px};") 
	f.write("\n")
	center_pts.append(j)
	elements.append(j)
	j=j+1
	
	#Right point
	f.write("Point("+str(j)+") = {center_x_"+str(i)+"+radius_"+str(i)+", center_y_"+str(i)+", height_"+str(i)+", volSize_px};")
	f.write("\n")
	right_pts.append(j)
	elements.append(j)
	j=j+1
	
	#Upper point
	f.write("Point("+str(j)+") = {center_x_"+str(i)+", center_y_"+str(i)+"+radius_"+str(i)+", height_"+str(i)+", volSize_px};")
	f.write("\n")
	upper_pts.append(j)
	elements.append(j)
	j=j+1
	
	#Left point
	f.write("Point("+str(j)+") = {center_x_"+str(i)+"-radius_"+str(i)+", center_y_"+str(i)+", height_"+str(i)+", volSize_px};")
	f.write("\n")
	left_pts.append(j)
	elements.append(j)
	j=j+1
	
	#Lower point
	f.write("Point("+str(j)+") = {center_x_"+str(i)+", center_y_"+str(i)+"-radius_"+str(i)+", height_"+str(i)+", volSize_px};")
	f.write("\n")
	lower_pts.append(j)
	elements.append(j)
	j=j+1
	
	#Empty vectors
	circle_temp=[] 
	
	#~~~~~~~~~~~
	#Arcs
	
	#Right top
	f.write("Circle("+str(j)+") = {"+str(right_pts[-1])+", "+ str(center_pts[-1]) +", "+ str(upper_pts[-1])+"};")
	f.write("\n")
	circle_temp.append(j)
	elements.append([right_pts[-1],center_pts[-1],upper_pts[-1]])
	j=j+1
	
	#Left top
	f.write("Circle("+str(j)+") = {"+str(upper_pts[-1])+", "+ str(center_pts[-1]) +", "+ str(left_pts[-1])+"};")
	f.write("\n")
	circle_temp.append(j)
	elements.append([upper_pts[-1],center_pts[-1],left_pts[-1]])
	j=j+1
	
	#Left bottom
	f.write("Circle("+str(j)+") = {"+str(left_pts[-1])+", "+ str(center_pts[-1]) +", "+ str(lower_pts[-1])+"};")
	f.write("\n")
	circle_temp.append(j)
	elements.append([left_pts[-1],center_pts[-1],lower_pts[-1]])
	j=j+1
	
	#Right bottom
	f.write("Circle("+str(j)+") = {"+str(lower_pts[-1])+", "+ str(center_pts[-1]) +", "+ str(right_pts[-1])+"};")
	f.write("\n")
	circle_temp.append(j)
	elements.append([lower_pts[-1],center_pts[-1],right_pts[-1]])
	j=j+1
	
	#Remebering circle for later loops and surfaces
	circles.append(circle_temp)

#------------------
#Create lines connecting circles

conn_line=[]
f.write("\n")
f.write("\n")
f.write("//Connecting lines ")
f.write("\n")

for i in range(1,len(files)):
	f.write("//Lines connecting stacks "+str(i-1)+"and " + str(i))
	f.write("\n")
	
	#~~~~~~~~~~~
	#Lines connecting circles
	conn_line_temp=[]
	
	f.write("Line("+str(j)+") = {"+str(right_pts[i])+","+str(right_pts[i-1])+"};")
	f.write("\n")
	conn_line_temp.append(j)
	elements.append([right_pts[i],right_pts[i-1]])
	j=j+1
	
	f.write("Line("+str(j)+") = {"+str(upper_pts[i])+","+str(upper_pts[i-1])+"};")
	f.write("\n")
	conn_line_temp.append(j)
	elements.append([upper_pts[i],upper_pts[i-1]])
	j=j+1	
		

	f.write("Line("+str(j)+") = {"+str(left_pts[i])+","+str(left_pts[i-1])+"};")
	f.write("\n")
	conn_line_temp.append(j)
	elements.append([left_pts[i],left_pts[i-1]])
	j=j+1
	
	
	f.write("Line("+str(j)+") = {"+str(lower_pts[i])+","+str(lower_pts[i-1])+"};")
	f.write("\n")
	conn_line_temp.append(j)
	elements.append([lower_pts[i],lower_pts[i-1]])
	j=j+1
	
	conn_line.append(conn_line_temp)
	
#------------------
#Create line loops and surfaces

f.write("\n")
f.write("\n")
f.write("//Surfaces ")
f.write("\n")
surfaces=[]
for i in range(1,len(files)):
	
	#Top right 
	f.write("//Top right connecting stacks "+str(i-1)+"and " + str(i))
	f.write("\n")
		
	conn_line[i-1][0],circles[i][0],conn_line[i-1][1],circles[i-1][0]=check_ends(elements,conn_line[i-1][0],circles[i][0],conn_line[i-1][1],circles[i-1][0])
	f.write("Line Loop("+str(j)+") = {"+str(conn_line[i-1][0])+"," +str(circles[i][0])+"," +str(conn_line[i-1][1])+"," +str(circles[i-1][0])+"};")
	j=j+1
	f.write("\n")
	f.write("Ruled Surface("+str(j)+") = {"+str(j-1)+"};")
	surfaces.append(j)
	j=j+1
	f.write("\n")
	
	#Top left
	f.write("//Top left connecting stacks "+str(i-1)+"and " + str(i))
	f.write("\n")
	
	conn_line[i-1][1],circles[i][1],conn_line[i-1][2],circles[i-1][1]=check_ends(elements,conn_line[i-1][1],circles[i][1],conn_line[i-1][2],circles[i-1][1])
	f.write("Line Loop("+str(j)+") = {"+str(conn_line[i-1][1])+"," +str(circles[i][1])+"," +str(conn_line[i-1][2])+"," +str(circles[i-1][1])+"};")
	j=j+1
	f.write("\n")
	f.write("Ruled Surface("+str(j)+") = {"+str(j-1)+"};")
	surfaces.append(j)
	j=j+1
	f.write("\n")
	
	#Bottom left
	f.write("//Top left connecting stacks "+str(i-1)+"and " + str(i))
	f.write("\n")
	
	conn_line[i-1][2],circles[i][2],conn_line[i-1][3],circles[i-1][2]=check_ends(elements,conn_line[i-1][2],circles[i][2],conn_line[i-1][3],circles[i-1][2])
	f.write("Line Loop("+str(j)+") = {"+str(conn_line[i-1][2])+"," +str(circles[i][2])+"," +str(conn_line[i-1][3])+"," +str(circles[i-1][2])+"};")
	j=j+1
	f.write("\n")
	f.write("Ruled Surface("+str(j)+") = {"+str(j-1)+"};")
	surfaces.append(j)
	j=j+1
	f.write("\n")
	
	#Bottom right 
	f.write("//Bottom right connecting stacks "+str(i-1)+"and " + str(i))
	f.write("\n")
	
	conn_line[i-1][3],circles[i][3],conn_line[i-1][0],circles[i-1][3]=check_ends(elements,conn_line[i-1][3],circles[i][3],conn_line[i-1][0],circles[i-1][3])
	f.write("Line Loop("+str(j)+") = {"+str(conn_line[i-1][3])+"," +str(circles[i][3])+"," +str(conn_line[i-1][0])+"," +str(circles[i-1][3])+"};")
	j=j+1
	f.write("\n")
	f.write("Ruled Surface("+str(j)+") = {"+str(j-1)+"};")
	surfaces.append(j)
	j=j+1
	f.write("\n")

#Top and bottom lid
f.write("//Top lid ")
f.write("\n")

f.write("Line Loop("+str(j)+") = {"+str(circles[0][0])+"," +str(circles[0][1])+"," +str(circles[0][2])+"," +str(circles[0][3])+"};")
j=j+1
f.write("\n")
f.write("Ruled Surface("+str(j)+") = {"+str(j-1)+"};")
surfaces.append(j)
j=j+1
f.write("\n")	

f.write("//Bottom lid ")
f.write("\n")

f.write("Line Loop("+str(j)+") = {"+str(circles[-1][0])+"," +str(circles[-1][1])+"," +str(circles[-1][2])+"," +str(circles[-1][3])+"};")
j=j+1
f.write("\n")
f.write("Ruled Surface("+str(j)+") = {"+str(j-1)+"};")
surfaces.append(j)
j=j+1
f.write("\n")	

#------------------
#Create surface loops

f.write("\n")
f.write("\n")
f.write("//Surface Loop ")
f.write("\n")
	
f.write("Surface Loop("+str(j)+") = {")	
for surf in surfaces:
	if surf !=surfaces[-1]:
		f.write(str(surf)+", ")	
	else:
		f.write(str(surf))	
f.write("};")
j=j+1
f.write("\n")

#------------------
#Create volumes

f.write("\n")
f.write("\n")
f.write("//Volumes ")
f.write("\n")

f.write("Volume("+str(j)+") = {"+str(j-1)+"};")

f.close()

#=========================================================================================================================
#Run Gmsh
#=========================================================================================================================

	
os.system("gmsh -3 " + folder+"/zStack.geo" )
		
	
	
	
	
	
	
	
	
	
	