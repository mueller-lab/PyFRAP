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
phi = CellVariable(name = "solution variable",
		mesh = mesh,
		value = 0.)

#Define term
D = 1.
eq = TransientTerm() == DiffusionTerm(coeff=D)

#Boundary and ICs
valueTopLeft = 0
valueBottomRight = 1

X, Y = mesh.faceCenters
facesTopLeft = ((mesh.facesLeft & (Y > L / 2))
		| (mesh.facesTop & (X < L / 2)))
facesBottomRight = ((mesh.facesRight & (Y < L / 2))
		| (mesh.facesBottom & (X > L / 2)))

phi.constrain(valueTopLeft, facesTopLeft)
phi.constrain(valueBottomRight, facesBottomRight)

#Plot
if __name__ == '__main__':
	viewer = Viewer(vars=phi, datamin=0., datamax=1.)
	viewer.plot()

#Time step
timeStepDuration = 10 * 0.9 * dx**2 / (2 * D)
steps = 10

for step in range(steps):
	eq.solve(var=phi,dt=timeStepDuration)

if __name__ == '__main__':
	viewer.plot()
