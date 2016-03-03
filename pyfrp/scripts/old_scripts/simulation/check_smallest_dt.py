#Script to run super simple 2D simulation taking from FiPy and check which is the smallest timestep to garanty no loss due to stepping

#Script to run super simple 2D simulation taking from FiPy

from fipy import *
from numpy import *
import matplotlib.pyplot as plt

#Some variables
cellSize = 0.05
radius = 1.
D = 1.
tend=1
	
#scales=[1,10,100,1000,10000,100000,1000000]
scales=[1,10,100]
#scales=[1,10]

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

tvec_all=[]
mphis_solve_all=[]
mphis_sweep_all=[]

for s in scales:
	print
	print "Starting with scale=", s
	
	#Variable
	phi_solve = CellVariable(name = "solution variable",mesh = mesh,value = 100.) 
	phi_sweep = CellVariable(name = "solution variable",mesh = mesh,value = 100.,hasOld=True) 
	
	#Define equation
	eq = TransientTerm() == DiffusionTerm(coeff=D)

	#Mesh centers
	x, y = mesh.cellCenters
	
	#Set IC
	phi_solve.setValue(0, where=(x < 0.3) & (x > -0.3) & (y < 0.3) & (y > -0.3))
	phi_sweep.setValue(0, where=(x < 0.3) & (x > -0.3) & (y < 0.3) & (y > -0.3))
		
	#Plot initial concentration
	#fig=plt.figure()
	#fig.show()
	#ax=fig.add_subplot(111)
	#ax.tricontourf(x,y,phi)
	#plt.draw()
	#raw_input()
	
	#Timestepping
	timeStepDuration = 10 * 0.9 * cellSize**2 / (2 * D)/s
	steps = int(ceil(tend/timeStepDuration))
	
	#Some output
	print "Going to perfrom ", steps, " timesteps."
	
	#Sweeping
	sweeps=10
	
	#Compute mean and sum of initial phi
	mphi_solve=mean(phi_solve.value)
	sphi_solve=sum(phi_solve.value)
	mphi_sweep=mean(phi_sweep.value)
	sphi_sweep=sum(phi_sweep.value)
	
	
	#Empty lists for bookkeeping
	mphis_solve=[]
	mphis_sweep=[]
	tvec=[]
	
	#Append initial values
	mphis_solve.append(mphi_solve)
	mphis_sweep.append(mphi_sweep)
	tvec.append(0)
	
	#Solve
	for step in range(steps):
		#Solve
		eq.solve(var=phi_solve,dt=timeStepDuration) 
	
		
		#print
		#print mean(phi.value)
		#print sum(phi.value)
		
		#Sweep
		for sweep in range(sweeps):
			phi_sweep.updateOld()
			res=eq.sweep(var=phi_sweep,dt=timeStepDuration)
			#print
			#print mean(phi.value)
			#print sum(phi.value)
		#print res	
		#raw_input()
		
		if mod(step,100)==0:
			print float(step)/float(steps)*100.
		
		#if mean(phi.value)-mphi!=0 or sum(phi.value)-sphi!=0:
			
			#print "Scale ", s, " too small"
			#print step,timeStepDuration , mphi,sphi
			#break
		
		#Update bookkeeping
		mphi_solve=mean(phi_solve.value)
		sphi_solve=sum(phi_solve.value)
		mphi_sweep=mean(phi_sweep.value)
		sphi_sweep=sum(phi_sweep.value)
	
		#Append to bookkeeping
		mphis_solve.append(mphi_solve)
		mphis_sweep.append(mphi_sweep)
		tvec.append(tvec[-1]+timeStepDuration)
	
	#Bookkeeping
	mphis_solve_all.append(mphis_solve)
	mphis_sweep_all.append(mphis_sweep)
	tvec_all.append(tvec)

#Plot mean concentration on normed mean concentration	
fig=plt.figure()
ax=fig.add_subplot(211)

for i in range(len(mphis_solve_all)):

	ax.plot(tvec_all[i],mphis_solve_all[i],c=((0,1,1/(i+1))),label="solve "+str(scales[i]))
	ax.plot(tvec_all[i],mphis_sweep_all[i],c=((1,1/(i+1),1)),label="sweep "+str(scales[i]))

plt.legend()

ax=fig.add_subplot(212)

for i in range(len(mphis_solve_all)):
	
	ax.plot(tvec_all[i],mphis_solve_all[i]/mphis_solve_all[i][0],c=((0,1,1/(i+1))),label="solve "+str(scales[i]))
	ax.plot(tvec_all[i],mphis_sweep_all[i]/mphis_sweep_all[i][0],c=((1,1/(i+1),1)),label="sweep "+str(scales[i]))
		 
plt.legend()
plt.show()
	
	


		
 
