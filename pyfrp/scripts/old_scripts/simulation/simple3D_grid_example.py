#Script to run super simple 3D diffusion simulation modified from FiPy on a regular 20x20x20 grid

#Import modules
from fipy import *

#Make grid
nx = 20
ny = nx
nz = ny
dx = 1.
dy = dx
dz = dx
L = dx * nx
mesh = Grid3D(dx=dx, dy=dy, dz=dz, nx=nx, ny=ny, nz=nz)

#Create variable
phi = CellVariable(name = "solution variable",mesh = mesh,value = 1.)

#Define term
D = 1.
eq = TransientTerm() == DiffusionTerm(coeff=D)

#Boundary and ICs
x,y,z=mesh.cellCenters
phi.setValue(0, where=(x < 5) & (x > 1) & (y < 5) & (y > 1))

sphis=[sum(phi.value)]
	
#Time step
timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 500

for step in range(steps):
	eq.solve(var=phi,dt=timeStepDuration)
	#print sum(phi.value)
	sphis.append(sum(phi.value))
	
print sphis[-1]/sphis[0]	
