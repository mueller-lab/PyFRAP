v
#=========================================================================================================================
#Automatically find connecting lines between stacks
#=========================================================================================================================
	
	
	
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
	f.write("height_"+str(i)+"="+str(i*stack_distance)+";")
	f.write("\n")

j=1

elements=[]
elements.append("Dummy")

stacks=[]
elmnt_to_coords=[]

#------------------
#Create splines per stack
for i in range(len(files)):
	
	f.write("//Spline of stack " + str(i))
	f.write("\n")
	
	mpoly=mpolys[i]
	
	stack=[]
	elmnt_to_coord=[]
	for k in range(len(mpoly)):
	
		#~~~~~~~~~~~
		#Points 
		f.write("Point("+str(j)+") = {"+str(mpoly[k,1])+", "+str(mpoly[k,0])+", height_"+str(i)+", volSize_px};") 
		f.write("\n")
		stack.append(j)
		elements.append(j)
		
		elmnt_to_coord.append([j,k])
		
		j=j+1
	
	stacks.append(stack)
	elmnt_to_coords.append(elmnt_to_coord)

#------------------
#Create stack connecting lines	
for i in range(len(stacks)-1):
	conn_lines=[]
	splines=[]
	lines=[]
	for k in range(shape(mpolys[i])[0]):
		ds=[]
		for k2 in range(shape(mpolys[i+1])[0]):
			ds.append(sqrt((mpolys[i][k,1]-mpolys[i+1][k2,1])**2+(mpolys[i][k,0]-mpolys[i+1][k2,0])**2)) 
			
		k2min=ds.index(min(ds))
		
		conn_line=[stacks[i][k],stacks[i+1][k2min]]
		conn_lines.append(conn_line)
		
		#~~~~~~~~~~~
		#Line
		f.write("Line("+str(j)+") = {"+str(conn_line[0])+","+str(conn_line[1])+"};") 
		f.write("\n")
		elements.append(j)
		lines.append(j)
		j=j+1
				
		
		if k>0:
			#~~~~~~~~~~~
			#Define two splines between two lines
			if conn_lines[-1][0]!=conn_lines[-2][0]:
				#Spline #1
				print "bla"
				
				
				f.write("Spline("+str(j)+") = {")		
			
				for l in range(conn_lines[-2][0],conn_lines[-1][0]):
					f.write(str(l)+", ")
				f.write(str(conn_lines[-1][0]+1)+"};") 
				f.write("\n")
				splines.append(j)
				elements.append(j)
				j=j+1	
			
			if conn_lines[-1][1]!=conn_lines[-2][1]:	
				#Spline #2
				f.write("Spline("+str(j)+") = {")		
				
				for l in range(conn_lines[-2][1],conn_lines[-1][1]):
					f.write(str(l)+", ")
				f.write(str(conn_lines[-1][1]+1)+"};") 
				f.write("\n")
				splines.append(j)
				elements.append(j)
				j=j+1
			
			
			print j-1,conn_lines[-1][0],conn_lines[-2][0],conn_lines[-1][1],conn_lines[-2][1]
				
			#~~~~~~~~~~~
			#Line Loop	
			#f.write("Line Loop("+str(j)+") = {"+str(lines[-2])+"," +str(splines[-1])+"," +str(lines[-1])+"," +str(splines)+"};")
		
		
		

	