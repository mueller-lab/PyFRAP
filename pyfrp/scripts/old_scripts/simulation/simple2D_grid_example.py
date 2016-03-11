#Script to run super simple 2D diffusion simulation taking from FiPy on a regular 20x20 grid

#Import modules
from fipy import *

#Make grid
nx = 20
ny = nx
dx = 1.
dy = dx
L = dx * nx
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

#Create variable
phi = CellVariable(name = "solution variable",mesh = mesh,value = 1.)

#Define term
D = 1.
eq = TransientTerm() == DiffusionTerm(coeff=D)

#Boundary and ICs
x,y=mesh.cellCenters
phi.setValue(0, where=(x < 5) & (x > 1) & (y < 5) & (y > 1))

print sum(phi.value)

#Plot
if __name__ == '__main__':
	viewer = Viewer(vars=phi, datamin=0., datamax=1.)
	viewer.plot()
	raw_input()
	
#Time step
timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 1000

for step in range(steps):
	eq.solve(var=phi,dt=timeStepDuration)
	print sum(phi.value)
if __name__ == '__main__':
	viewer.plot()
	raw_input()