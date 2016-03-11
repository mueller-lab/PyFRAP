#Script to run super simple 3D simulation, and test the influence of mesh sizes onto the gibbs effect

#Import some modules
from fipy import *
from numpy import *
import matplotlib.pyplot as plt
import matplotlib as mpl

import sys

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

def profile_plots(ax,x,line,line_ind,prof_every,bins=False,nbins=100):
	
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
			
			#Convert to np array
			x_sort=asarray(x_sort)
			phi_sort=asarray(phi_sort)
			
				
			#Plot
			c = mpl.cm.jet(float(k)/(float(kmax)),1)
			
			if bins:
				bins_x=linspace(min(x_sort)+0.01*min(x_sort),max(x_sort)+0.01*max(x_sort),nbins)
				#digitized = digitize(phi_sort, bins_x)
				#bin_means = [phi_sort[digitized == j].mean() for j in range(1, len(bins_x))]
				#bins_x = [(bins_x[j+1] + bins_x[j])/2. for j in range(len(bins_x)-1)]
				#bins_x=bins_x[:-1]
				bins_mid,hist_phi,bin_phi=simple_hist(bins_x,x_sort,phi_sort)
				
				print len(bins_mid), len(bin_phi)
				ax.plot(bins_mid,bin_phi,c=c,label="step = "+str(i))
			else:
				ax.plot(x_sort,phi_sort,c=c,label="step = "+str(i))
		
			k=k+1
	
	plt.legend(bbox_to_anchor=(1.3, 1))
	
	return ax

def comp_error(vec):
	return min(vec)/max(vec)

def get_all_errors(squ,out,sl,mphis):
	err_squ=comp_error(squ)
	err_out=comp_error(out)
	err_sl=comp_error(sl)
	err_mphis=comp_error(mphis)
	
	return err_squ, err_out, err_sl, err_mphis
	
def simple_hist(bins,data_x,data_y,is_sorted=False):
	hist_y=zeros(len(bins)-1)
	bin_y=zeros(len(bins)-1)
	
	for i in range(len(bins)-1):
		for j in range(len(data_x)):
			if i==0:
				if bins[i]<= data_x[j] and data_x[j]<=bins[i+1]:
					bin_y[i]=bin_y[i]+data_y[j]
					hist_y[i]=hist_y[i]+1
			else:
				if bins[i]< data_x[j] and data_x[j]<=bins[i+1]:
					bin_y[i]=bin_y[i]+data_y[j]
					hist_y[i]=hist_y[i]+1
	
	bins_mid = [(bins[j+1] + bins[j])/2. for j in range(len(bins)-1)]				
	
	return bins_mid, hist_y, bin_y

def compl_inds(l1,l2):
	
	l=[]
	for i in l1:
		if i not in l2:
			l.append(i)
	
	return l

			
def threeD_sim(cellSize=0.05,radius=1.,height=1.,D=1.,val_in=0.,val_out=1.,prof_every=100,steps=500,sweeps=3,method='sweep',liveplot=False,avg='both',norm=True,dt_scale=0.01,debug=False,addNoise=False,sigmaNoise=0.1,res_req=inf,final_plot=True,mesh=None):

	#Mesh
	if mesh==None:
		mesh = Gmsh3D(''' cellSize = %(cellSize)g; 
			radius = %(radius)g;
			height = %(height)g; 
			Point(1) = {0, 0, 0, cellSize}; 
			Point(2) = {-radius, 0, 0, cellSize};
			Point(3) = {0, radius, 0, cellSize};
			Point(4) = {radius, 0, 0, cellSize};
			Point(5) = {0, -radius, 0, cellSize};
			Circle(6) = {2, 1, 3}; 
			Circle(7) = {3, 1, 4}; 
			Circle(8) = {4, 1, 5}; 
			Circle(9) = {5, 1, 2};
			Line Loop(10) = {6,7,8,9};
			Plane Surface(11) = {10};
			Extrude {0,0,height} {
				Surface{11};
			}''' % locals()) 
	
	#Variable
	if method=='sweep':
		phi = CellVariable(name = "solution variable",mesh = mesh,value = val_out,hasOld=True) 
	elif method=='solve':
		phi = CellVariable(name = "solution variable",mesh = mesh,value = 1.) 

	#Define equation
	eq = TransientTerm() == DiffusionTerm(coeff=D)

	#Mesh centers
	x, y , z = mesh.cellCenters

	#Get getCellVolumes
	cvs=mesh.getCellVolumes()

	#Set IC
	phi.setValue(val_in, where=(x < 0.3) & (x > -0.3) & (y < 0.3) & (y > -0.3))
	
	#Add add noise 
	if addNoise:
		n=random.normal(0,sigmaNoise,len(x))
		phi.value=phi.value+n
		
	#Backup IC values
	phi_vals_IC=phi.value.copy()

	#Get indices of all regions
	slice_ind=where((z<0.7*height) & (z > 0.3*height))[0]
	squ_ind=where((x<0.3) & (x > - 0.3) & (y < 0.3) & (y > -0.3) & (z<0.7*height) & (z > 0.3*height))[0]
	left_ind=where(x<0)[0]
	right_ind=where(x>0)[0]
	out_ind=compl_inds(slice_ind,squ_ind)
	line_ind=where((y < 0.1) & (y > -0.1))[0]
	all_ind=range(len(x))
	
	#Grab x,y coordinates of slice
	x_slice=x[slice_ind]
	y_slice=y[slice_ind]
	phi_slice=phi.value[slice_ind]
	
	#Timestepping
	#timeStepDuration = 10 * 0.9 * cellSize**2 / (2 * D)*dt_scale
	timeStepDuration =  0.9  / (2 * D)*dt_scale
	
	
	#Track mean/sum of phi
	mphis=[]
	sphis=[]
	tvec=[]
	mphis_abs=[]
	sphis_abs=[]
	
	mphis.append(get_conc(phi,cvs,all_ind))
	sphis.append(sum(phi.value*cvs))
	mphis_abs.append(get_abs(phi,all_ind))
	sphis_abs.append(sum(phi.value))
	tvec.append(0)
	
	#Keep track of phi itself
	phis=[phi.value]
	
	#Keep track of sweeps and res
	ress=[]
	ress0=[]
	sweepss=[]
	
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
	mySolver = LinearLUSolver(iterations=1234511111, tolerance=5e-20)
	#mySolver = LinearPCGSolver(iterations=1234, tolerance=5e-6)
	
	#Print some debugging messages
	if debug:
		print "============================================================================================================================================"
		print "====Mesh===="
		print "ncells = ", len(cvs)
		print "min(cvs) = ", min(cvs)
		print "max(cvs) = ", max(cvs)
		print "mean(cvs) = ", mean(cvs)
		print min(x)," < x < ", max(x)
		print min(y)," < y < ", max(y)
		print min(z)," < z < ", max(z)
		print
		
		print "ncells[squ_ind] = ", len(cvs[squ_ind])
		print "min(cvs[squ_ind]) = ", min(cvs[squ_ind])
		print "max(cvs[squ_ind]) = ", max(cvs[squ_ind])
		print "mean(cvs[squ_ind]) = ", mean(cvs[squ_ind])
		print min(x[squ_ind])," < x < ", max(x[squ_ind])
		print min(y[squ_ind])," < y < ", max(y[squ_ind])
		print min(z[squ_ind])," < z < ", max(z[squ_ind])
		print
		
		print "ncells[out_ind] = ", len(cvs[out_ind])
		print "min(cvs[out_ind]) = ", min(cvs[out_ind])
		print "max(cvs[out_ind]) = ", max(cvs[out_ind])
		print "mean(cvs[out_ind]) = ", mean(cvs[out_ind])
		print min(x[out_ind])," < x < ", max(x[out_ind])
		print min(y[out_ind])," < y < ", max(y[out_ind])
		print min(z[out_ind])," < z < ", max(z[out_ind])
		print
		
		print "ncells[slice_ind] = ", len(cvs[slice_ind])
		print "min(cvs[slice_ind]) = ", min(cvs[slice_ind])
		print "max(cvs[slice_ind]) = ", max(cvs[slice_ind])
		print "mean(cvs[slice_ind]) = ", mean(cvs[slice_ind])
		print min(x[slice_ind])," < x < ", max(x[slice_ind])
		print min(y[slice_ind])," < y < ", max(y[slice_ind])
		print min(z[slice_ind])," < z < ", max(z[slice_ind])
		print
		
		print "ncells[right_ind] = ", len(cvs[right_ind])
		print "min(cvs[right_ind]) = ", min(cvs[right_ind])
		print "max(cvs[right_ind]) = ", max(cvs[right_ind])
		print "mean(cvs[right_ind]) = ", mean(cvs[right_ind])
		print min(x[right_ind])," < x < ", max(x[right_ind])
		print min(y[right_ind])," < y < ", max(y[right_ind])
		print min(z[right_ind])," < z < ", max(z[right_ind])
		print
		
		print "ncells[left_ind] = ", len(cvs[left_ind])
		print "min(cvs[left_ind]) = ", min(cvs[left_ind])
		print "max(cvs[left_ind]) = ", max(cvs[left_ind])
		print "mean(cvs[left_ind]) = ", mean(cvs[left_ind])
		print min(x[left_ind])," < x < ", max(x[left_ind])
		print min(y[left_ind])," < y < ", max(y[left_ind])
		print min(z[left_ind])," < z < ", max(z[left_ind])
		print
		
	
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
		axes,backgrounds,p=fast_cont(x_slice,y_slice,phi_slice,figlive,axes,backgrounds,clevels=clevels,cmap=cmap)

		#Draw colorbar, save backgrounds again and redraw
		#plt.colorbar(p)
		#backgrounds = [figlive.canvas.copy_from_bbox(ax.bbox) for ax in axes]
		#figlive.canvas.draw()
		raw_input()
	
	#Solve
	for step in range(steps):
		
		if method=="sweep":
			sweep=0
			res=inf
			phi.updateOld()
			eq.solve(var=phi,dt=timeStepDuration,solver=mySolver)
			while res>res_req or sweep < sweeps:
				
				res=eq.sweep(var=phi,dt=timeStepDuration,solver=mySolver)
				
				if sweep==0:
					ress0.append(res)
				
				sweep=sweep+1
				
				
				
				#print res,sweep
			#raw_input()		
			sweepss.append(sweep)
			ress.append(res)
		
		elif method=="solve":		
			eq.solve(var=phi,dt=timeStepDuration,solver=mySolver) 
		
		curr_perc=int(float(step)/(steps)*100)
		sys.stdout.write("\r%d%%" %curr_perc)  
		sys.stdout.flush()
		
		#Track mean/sum of phi
		mphis.append(get_conc(phi,cvs,all_ind))
		sphis.append(sum(phi.value*cvs))
		mphis_abs.append(get_abs(phi,all_ind))
		sphis_abs.append(sum(phi.value))
		tvec.append(tvec[-1]+timeStepDuration)
		
		#Keep track of phi
		phis.append(phi.value)
		
		#Grab x,y coordinates of slice
		x_slice=x[slice_ind]
		y_slice=y[slice_ind]
		phi_slice=phi.value[slice_ind]
		
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
	if final_plot:
		fig=plt.figure()

		#Normed sum(phi)/mean(phi)
		ax=fig.add_subplot(231)
		ax.set_title("Conservation Check")
		ax.plot(tvec,sphis/sphis[0],'m',label='sum(phi*cvs) normed')
		ax.plot(tvec,mphis/mphis[0],'c',label='mean(phi*cvs) normed')
		ax.plot(tvec,sphis_abs/sphis_abs[0],'r',label='sum(phi) normed')
		ax.plot(tvec,mphis_abs/mphis_abs[0],'b',label='mean(phi) normed')
		
		ax.set_ylim([0,2])
		plt.legend()

		#Region specific concentrations absolute
		if avg in ['abs','both']:
			ax=fig.add_subplot(232)
			ax.set_title("Phi[regions] (concentrations)")
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
				
			plt.legend(loc=4)
		
		#Region specific concentrations concentration
		if avg in ['conc','both']:
			ax=fig.add_subplot(233)
			ax.set_title("Phi*cvs[regions] (absolute value)")
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
			
			plt.legend(loc=4)
		
		#Final contour
		ax=fig.add_subplot(235)
		ax.set_title("Phi(t=tend)")
		clevels_IC=linspace(val_in,val_out,50)
		a=ax.tricontourf(x_slice,y_slice,phi_slice,levels=clevels_IC)
		plt.colorbar(a)

		#Profile
		ax.set_title("Phi(-0.1<y<0.1)")
		ax=fig.add_subplot(236)
		ax=profile_plots(ax,x,line,line_ind,prof_every)
				
		#cvs contour
		ax=fig.add_subplot(234)
		ax.set_title("CVS[slice_ind]")
		ax.tricontourf(x_slice,y_slice,cvs[slice_ind])
		plt.colorbar(a)		
		
		#Figure 2
		#Plot distribution of phi over different cvs
		fig2=plt.figure()
		ax=fig2.add_subplot(231)
		profile_plots(ax,cvs,phis,range(len(cvs)),prof_every,bins=True,nbins=100)
		
		#Plot location of outer indices
		ax=fig2.add_subplot(232)
		out_reg=zeros(shape(x))
		out_reg[out_ind]=1
		ax.scatter(x[out_ind],y[out_ind])
		
		#Plot location of squ indices
		ax=fig2.add_subplot(233)
		squ_reg=zeros(shape(x))
		squ_reg[squ_ind]=1
		ax.scatter(x[squ_ind],y[squ_ind])
		
		if method=="sweep":
			
			#Plot number of sweeps taken to accomplish end criteria
			ax=fig2.add_subplot(234)
			ax.plot(tvec[1:],sweepss,'b',label='sweeps')
			sweeps_max=sweeps*ones(shape(tvec[1:]))
			ax.plot(tvec[1:],sweeps_max,'r',label="sweeps max")	
			ax.set_title("Sweeps required")
			ax.set_ylim([min([min(sweepss),sweeps])*0.7,max([max(sweepss),sweeps])*1.1])
			plt.legend()
			
			#Plot starting and final res
			ax=fig2.add_subplot(235)
			ax.plot(tvec[1:],ress,'b',label="res")
			ax.plot(tvec[1:],ress0,'g',label="res0")
			
			if not isinf(res_req):
				
				reqs=res_req*ones(shape(tvec[1:]))
				ax.plot(tvec[1:],reqs,'r',label="res_req")
				ax.set_ylim([min([min(ress),res_req,min(ress0)])*0.9,max([max(ress),res_req,max(ress0)])*1.1])
				
			ax.set_title("Final res")
			plt.legend()
		
	#Print final Conservation statistics
	if debug:
		print "====Numerical errors===="
		print "sphis[0]/sphis[-1]:", sphis[0]/sphis[-1]
		print "sphis_abs[0]/sphis_abs[-1]",sphis_abs[0]/sphis_abs[-1]
		print "min(sphis)/max(sphis)",min(sphis)/max(sphis)
		print "min(sphis_abs)/max(sphis_abs)",min(sphis_abs)/max(sphis_abs)
		print "Timepoint of smallest value:", sphis_abs.index(min(sphis_abs))
		print "Timepoint of largest value:", sphis_abs.index(max(sphis_abs))
		print 
		
		print "mphis[0]/mphis[-1]:", mphis[0]/mphis[-1]
		print "mphis_abs[0]/mphis_abs[-1]",mphis_abs[0]/mphis_abs[-1]
		print "min(mphis)/max(mphis)",min(mphis)/max(mphis)
		print "min(mphis_abs)/max(mphis_abs)",min(mphis_abs)/max(mphis_abs)
		print "Timepoint of smallest value:", mphis_abs.index(min(mphis_abs))
		print "Timepoint of largest value:", mphis_abs.index(max(mphis_abs))
		print 
		
		print "sl_conc[0]/sl_conc[-1]:", sl_conc[0]/sl_conc[-1]
		print "sl_abs[0]/sl_abs[-1]", sl_abs[0]/sl_abs[-1]
		print "min(sl)/max(sl)",min(sl_conc)/max(sl_conc)
		print "min(sl_abs)/max(sl_abs)",min(sl_abs)/max(sl_abs)
		print "Timepoint of smallest value:", sl_abs.index(min(sl_abs))
		print "Timepoint of largest value:", sl_abs.index(max(sl_abs))
		print 
		
		print "out_conc[0]/out_conc[-1]:", out_conc[0]/out_conc[-1]
		print "out_abs[0]/out_abs[-1]", out_abs[0]/out_abs[-1]
		print "min(out)/max(out)",min(out_conc)/max(out_conc)
		print "min(out_abs)/max(out_abs)",min(out_abs)/max(out_abs)
		print "Timepoint of smallest value:", out_abs.index(min(out_abs))
		print "Timepoint of largest value:", out_abs.index(max(out_abs))
		print 
		
		print "squ_conc[0]/squ_conc[-1]:", squ_conc[0]/squ_conc[-1]
		print "squ_abs[0]/squ_abs[-1]", squ_abs[0]/squ_abs[-1]
		print "min(squ)/max(squ)",min(squ_conc)/max(squ_conc)
		print "min(squ_abs)/max(squ_abs)",min(squ_abs)/max(squ_abs)
		print "Timepoint of smallest value:", squ_abs.index(min(squ_abs))
		print "Timepoint of largest value:", squ_abs.index(max(squ_abs))
		print 
		
		print "right_conc[0]/right_conc[-1]:", right_conc[0]/right_conc[-1]
		print "right_abs[0]/right_abs[-1]", right_abs[0]/right_abs[-1]
		print "min(right)/max(right)",min(right_conc)/max(right_conc)
		print "min(right_abs)/max(right_abs)",min(right_abs)/max(right_abs)
		print "Timepoint of smallest value:", right_abs.index(min(right_abs))
		print "Timepoint of largest value:", right_abs.index(max(right_abs))
		print 
		
		print "left_conc[0]/left_conc[-1]:", left_conc[0]/left_conc[-1]
		print "left_abs[0]/left_abs[-1]", left_abs[0]/left_abs[-1]
		print "min(left)/max(left)",min(left_conc)/max(left_conc)
		print "min(left_abs)/max(left_abs)",min(left_abs)/max(left_abs)
		print "Timepoint of smallest value:", left_abs.index(min(left_abs))
		print "Timepoint of largest value:", left_abs.index(max(left_abs))
		print 
		
		#raw_input("Done, Press Enter")
	
	
	if final_plot:
		plt.show()
	
	return tvec,squ_abs,out_abs,sl_abs,mphis_abs,squ_conc,out_conc,sl_conc,mphis,mesh


#Main script
steps=250
dt_scale=0.001
debug=True

#mesh_sizes=arange(0.08,0.1,0.5,0.05)
#mesh_sizes=[0.07,0.08,0.09,0.1,0.25,0.5]
mesh_sizes=[0.15]

fig=plt.figure()
ax_squ=fig.add_subplot(221)
ax_out=fig.add_subplot(222)
ax_sl=fig.add_subplot(223)
ax_m=fig.add_subplot(224)

fig2=plt.figure()
ax_squ2=fig2.add_subplot(221)
ax_out2=fig2.add_subplot(222)
ax_sl2=fig2.add_subplot(223)
ax_m2=fig2.add_subplot(224)

#Loop through different mesh sizes
kmax=len(mesh_sizes)

#Some bookkeeping
errs_squ=[]
errs_out=[]
errs_sl=[]
errs_left=[]
errs_right=[]
errs_mphis=[]

errs_squ_abs=[]
errs_out_abs=[]
errs_sl_abs=[]
errs_left_abs=[]
errs_right_abs=[]
errs_mphis_abs=[]

ncells=[]

#Run with res_req=0.5
for k,cellSize in enumerate(mesh_sizes):

	#Run with solve
	tvec,squ_abs,out_abs,sl_abs,mphis_abs,squ_conc,out_conc,sl_conc,mphis_conc,mesh=threeD_sim(method='solve',liveplot=False,steps=steps,prof_every=50,norm=False,cellSize=cellSize,sweeps=5,dt_scale=dt_scale,debug=debug,addNoise=False,res_req=0.05,final_plot=True)
		
	print "Done with cellSize = ", cellSize
	
	#Plot latest timeseries
	c = mpl.cm.jet(float(k)/(float(kmax)),1)
	
	ax_squ.plot(tvec,squ_abs,c=c,label=str(cellSize))
	ax_out.plot(tvec,out_abs,c=c,label=str(cellSize))
	ax_sl.plot(tvec,sl_abs,c=c,label=str(cellSize))
	ax_m.plot(tvec,mphis_abs,c=c,label=str(cellSize))

	ax_squ2.plot(tvec,squ_conc,c=c,label=str(cellSize))
	ax_out2.plot(tvec,out_conc,c=c,label=str(cellSize))
	ax_sl2.plot(tvec,sl_conc,c=c,label=str(cellSize))
	ax_m2.plot(tvec,mphis_conc,c=c,label=str(cellSize))
	
	#Compute Errors
	err_squ, err_out, err_sl, err_mphis = get_all_errors(squ_conc,out_conc,sl_conc,mphis_conc)
	errs_squ.append(err_squ)
	errs_out.append(err_out)
	errs_sl.append(err_sl)
	errs_mphis.append(err_mphis)
	
	err_squ, err_out, err_sl, err_mphis = get_all_errors(squ_abs,out_abs,sl_abs,mphis_abs)
	errs_squ_abs.append(err_squ)
	errs_out_abs.append(err_out)
	errs_sl_abs.append(err_sl)
	errs_mphis_abs.append(err_mphis)
	
	#Keep track of mesh size
	ncells.append(len(mesh.x))
	
ax_out.set_title('out')
ax_sl.set_title('slice')
ax_m.set_title('all')
ax_squ.set_title('squ')

plt.legend()

ax_out2.set_title('out')
ax_sl2.set_title('slice')
ax_m2.set_title('all')	
ax_squ2.set_title('squ')

plt.legend()

fig3=plt.figure()
ax_err_squ=fig3.add_subplot(321)
ax_err_out=fig3.add_subplot(322)
ax_err_sl=fig3.add_subplot(323)
ax_err_mphis=fig3.add_subplot(324)
ax_ncells=fig3.add_subplot(325)

fig4=plt.figure()
ax_err_squ_abs=fig4.add_subplot(321)
ax_err_out_abs=fig4.add_subplot(322)
ax_err_sl_abs=fig4.add_subplot(323)
ax_err_mphis_abs=fig4.add_subplot(324)
ax_ncells2=fig4.add_subplot(325)

c='r'
	
ax_err_squ.plot(mesh_sizes,errs_squ,c=c,label=str(cellSize))
ax_err_out.plot(mesh_sizes,errs_out,c=c,label=str(cellSize))
ax_err_sl.plot(mesh_sizes,errs_sl,c=c,label=str(cellSize))
ax_err_mphis.plot(mesh_sizes,errs_mphis,c=c,label=str(cellSize))
ax_ncells.plot(mesh_sizes,ncells,c=c,label=str(cellSize))

ax_err_out.set_title('out err conc')
ax_err_sl.set_title('slice err conc')
ax_err_mphis.set_title('all err conc')
ax_err_squ.set_title('squ err conc')
ax_ncells.set_title('ncells')

ax_err_squ_abs.plot(mesh_sizes,errs_squ_abs,c=c,label=str(cellSize))
ax_err_out_abs.plot(mesh_sizes,errs_out_abs,c=c,label=str(cellSize))
ax_err_sl_abs.plot(mesh_sizes,errs_sl_abs,c=c,label=str(cellSize))
ax_err_mphis_abs.plot(mesh_sizes,errs_mphis_abs,c=c,label=str(cellSize))
ax_ncells2.plot(mesh_sizes,ncells,c=c,label=str(cellSize))

ax_err_out_abs.set_title('out err abs')
ax_err_sl_abs.set_title('slice err abs')
ax_err_mphis_abs.set_title('all err abs')
ax_err_squ_abs.set_title('squ  err abs')
ax_ncells2.set_title('ncells')



plt.show()