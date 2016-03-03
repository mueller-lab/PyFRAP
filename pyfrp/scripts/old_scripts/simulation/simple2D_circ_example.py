#Script to run super simple 2D simulation taking from FiPy

#Import some modules
from fipy import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib as mpl

def get_conc(var,cvs,ind):
	return sum(var.value[ind]*cvs[ind])/sum(cvs[ind])

def get_abs(var,ind):
	return mean(var.value[ind])

def fast_cont(x,y,var,fig,axes,backgrounds,clevels=None,cmap=mpl.cm.jet):
	
	items = enumerate(zip(axes, backgrounds), start=1)
	for j, (ax, background) in items:
		
		fig.canvas.restore_region(background)
		
		p = ax.tricontourf(x,y,var,levels=clevels,cmap=cmap)
		for patch in p.collections:
			ax.draw_artist(patch)
		fig.canvas.blit(ax.bbox)
		
	return axes,backgrounds,p

def profile_plots(ax,x,line,line_ind,prof_every):
	
	k=0
	kmax=len(line)/prof_every+1
	for i in range(len(line)):
		if mod(i,prof_every)==0:

			#Sort with respect to x
			xline=x[line_ind]
			s=sorted(zip(xline, line[i]), key=lambda xline:xline[0])
			x_sort=[]
			phi_sort=[]
			
			for j in range(shape(s)[0]):
				x_sort.append(s[j][0])
				phi_sort.append(s[j][1])
				
			#Plot
			c = mpl.cm.jet(float(k)/(float(kmax)),1)
			ax.plot(x_sort,phi_sort,c=c,label="step = "+str(i))
		
			k=k+1
	
	plt.legend()
	
	return ax
			
def twoD_sim(cellSize=0.5,radius=1.,D=1.,val_in=0.,val_out=1.,prof_every=100,steps=500,sweeps=3,method='sweep',liveplot=False,avg='both',norm=True):

	#Mesh
	mesh = Gmsh2D('''
		cellSize = %(cellSize)g;
		radius = %(radius)g;
		Point(1) = {0, 0, 0, cellSize};
		Point(2) = {-radius, 0, 0, cellSize};
		Point(3) = {0, radius, 0, cellSize};
		Point(4) = {radius, 0, 0, cellSize};
		Point(5) = {0, -radius, 0, cellSize};
		Circle(6) = {2, 1, 3};
		Circle(7) = {3, 1, 4};
		Circle(8) = {4, 1, 5};
		Circle(9) = {5, 1, 2};
		Line Loop(10) = {6, 7, 8, 9};
		Plane Surface(11) = {10};
		''' % locals()) 

	#Variable
	if method=='sweep':
		phi = CellVariable(name = "solution variable",mesh = mesh,value = val_out,hasOld=True) 
	elif method=='solve':
		phi = CellVariable(name = "solution variable",mesh = mesh,value = 1.) 

	#Define equation
	eq = TransientTerm() == DiffusionTerm(coeff=D)

	#Mesh centers
	x, y = mesh.cellCenters

	#Get getCellVolumes
	cvs=mesh.getCellVolumes()

	#Set IC
	phi.setValue(val_in, where=(x < 0.3) & (x > -0.3) & (y < 0.3) & (y > -0.3))

	#Backup IC values
	phi_vals_IC=phi.value.copy()

	#Get indices of all regions
	slice_ind=range(len(x))
	squ_ind=where((x<0.3) & (x > - 0.3) & (y < 0.3) & (y > -0.3))[0]
	left_ind=where(x<0)[0]
	right_ind=where(x>0)[0]
	out_ind=delete(slice_ind,squ_ind)
	line_ind=where((y < 0.1) & (y > -0.1))[0]

	#Timestepping
	timeStepDuration = 10 * 0.9 * cellSize**2 / (2 * D)*0.001

	#Compute mean/sum of phi
	mphi=mean(phi.value*cvs)
	sphi=sum(phi.value*cvs)

	#Track mean/sum of phi
	mphis=[]
	sphis=[]
	tvec=[]

	mphis.append(mphi)
	sphis.append(sphi)
	tvec.append(0)

	#Track region specific concentrations
	
	if avg in ['conc','both']:
		squ_conc=[]
		out_conc=[]
		sl_conc=[]
		left_conc=[]
		right_conc=[]
		
		squ_conc.append(get_conc(phi,cvs,squ_ind))
		out_conc.append(get_conc(phi,cvs,out_ind))
		sl_conc.append(get_conc(phi,cvs,slice_ind))
		left_conc.append(get_conc(phi,cvs,right_ind))
		right_conc.append(get_conc(phi,cvs,left_ind))
	
	if avg in ['abs','both']:
		squ_abs=[]
		out_abs=[]
		sl_abs=[]
		left_abs=[]
		right_abs=[]
		
		squ_abs.append(get_abs(phi,squ_ind))
		out_abs.append(get_abs(phi,out_ind))
		sl_abs.append(get_abs(phi,slice_ind))
		left_abs.append(get_abs(phi,right_ind))
		right_abs.append(get_abs(phi,left_ind))
	
	#Line profile
	line=[]
	line.append(phi.value[line_ind])

	#Create solver
	mySolver = LinearLUSolver(iterations=12345111, tolerance=5e-15)
	#mySolver = LinearPCGSolver(iterations=1234, tolerance=5e-6)

	#Live plotting
	if liveplot:
	
		figlive=plt.figure()
		figlive.show()
		ax=figlive.add_subplot(111)
		axes=[ax]

		#Make colorbar
		clevels=linspace(phi.value.min(),phi.value.max(),25)
		cmap = mpl.cm.jet

		# Let's capture the background of the figure
		backgrounds = [figlive.canvas.copy_from_bbox(ax.bbox) for ax in axes]

		# We need to draw the canvas before we start animating...
		figlive.canvas.draw()

		#Fast plot contour
		axes,backgrounds,p=fast_cont(x,y,phi.value,figlive,axes,backgrounds,clevels=clevels,cmap=cmap)

		#Draw colorbar, save backgrounds again and redraw
		plt.colorbar(p)
		backgrounds = [figlive.canvas.copy_from_bbox(ax.bbox) for ax in axes]
		figlive.canvas.draw()

	#Solve
	for step in range(steps):
		
		if method=="sweep":
			for sweep in range(sweeps):
				phi.updateOld()
				res=eq.sweep(var=phi,dt=timeStepDuration,solver=mySolver) 
		elif method=="solve":		
			eq.solve(var=phi,dt=timeStepDuration,solver=mySolver) 
		
		if mod(step,prof_every)==0:
			print step
		
		#Average 
		mphi=mean(phi.value*cvs)
		sphi=sum(phi.value*cvs)
		
		#Track mean/sum of phi
		mphis.append(mphi)
		sphis.append(sphi)
		tvec.append(tvec[-1]+timeStepDuration)
		
		#Track region specific concentrations
		if avg in ['conc','both']:
			squ_conc.append(get_conc(phi,cvs,squ_ind))
			out_conc.append(get_conc(phi,cvs,out_ind))
			sl_conc.append(get_conc(phi,cvs,slice_ind))
			left_conc.append(get_conc(phi,cvs,right_ind))
			right_conc.append(get_conc(phi,cvs,left_ind))
	
		if avg in ['abs','both']:
			squ_abs.append(get_abs(phi,squ_ind))
			out_abs.append(get_abs(phi,out_ind))
			sl_abs.append(get_abs(phi,slice_ind))
			left_abs.append(get_abs(phi,right_ind))
			right_abs.append(get_abs(phi,left_ind))
		
		#Fast plot contour
		if liveplot:
			axes,backgrounds,p=fast_cont(x,y,phi.value,figlive,axes,backgrounds,clevels=clevels)
		
		#Profile
		line.append(phi.value[line_ind])
		
	#Final plot	
	fig=plt.figure()

	#Normed sum(phi)/mean(phi)
	ax=fig.add_subplot(231)
	ax.plot(tvec,sphis/sphis[0],'m',label='sum(phi) normed')
	ax.plot(tvec,mphis/mphis[0],'c',label='mean(phi) normed')
	ax.set_title('')
	plt.legend()

	#Region specific concentrations absolute
	if avg in ['abs','both']:
		ax=fig.add_subplot(232)
		if norm:
			ax.plot(tvec,squ_abs/sl_abs[0],'b',label='squ')
			ax.plot(tvec,out_abs/sl_abs[0],'r',label='out')
			ax.plot(tvec,sl_abs/sl_abs[0],'g',label='slice')
			ax.plot(tvec,right_abs/sl_abs[0],'m',label='right')
			ax.plot(tvec,left_abs/sl_abs[0],'c',label='left')
		else:
			ax.plot(tvec,squ_abs,'b',label='squ')
			ax.plot(tvec,out_abs,'r',label='out')
			ax.plot(tvec,sl_abs,'g',label='slice')
			ax.plot(tvec,right_abs,'m',label='right')
			ax.plot(tvec,left_abs,'c',label='left')
			
		plt.legend()
	
	#Region specific concentrations concentration
	if avg in ['conc','both']:
		ax=fig.add_subplot(233)
		if norm:
			ax.plot(tvec,squ_conc/sl_conc[0],'b',label='squ')
			ax.plot(tvec,out_conc/sl_conc[0],'r',label='out')
			ax.plot(tvec,sl_conc/sl_conc[0],'g',label='slice')
			ax.plot(tvec,right_conc/sl_conc[0],'m',label='right')
			ax.plot(tvec,left_conc/sl_conc[0],'c',label='left')
		else:
			ax.plot(tvec,squ_conc,'b',label='squ')
			ax.plot(tvec,out_conc,'r',label='out')
			ax.plot(tvec,sl_conc,'g',label='slice')
			ax.plot(tvec,right_conc,'m',label='right')
			ax.plot(tvec,left_conc,'c',label='left')
		
		plt.legend()
	
	#Final contour
	ax=fig.add_subplot(235)
	clevels_IC=linspace(val_in,val_out,50)
	a=ax.tricontourf(x,y,phi.value,levels=clevels_IC)
	plt.colorbar(a)

	#Profile
	ax=fig.add_subplot(236)
	ax=profile_plots(ax,x,line,line_ind,prof_every)
			
	#cvs contour
	ax=fig.add_subplot(234)
	a=ax.tricontourf(x,y,cvs)
	plt.colorbar(a)		
			
	plt.show()

twoD_sim(method='solve',liveplot=True,steps=750,prof_every=50,norm=False,cellSize=0.25)

