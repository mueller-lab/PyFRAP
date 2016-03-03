###Script to debug embryo/mesh class on cluster

#Import module
import pyfrp_embryo

#Output meshfile
fnMsh='meshfiles/cylinder.msh'

#Geometrical parameters
center=[256,256]
cylinderHeight=90.3332804037
imagingRadius=294.103867936	
	
#Create embryo	
emb=pyfrp_embryo.embryo("TestedOnCluster")	

#Create geometry
emb.setGeometry2Cylinder(center,imagingRadius,cylinderHeight)
	
#Create Simulation
sim=emb.newSimulation()	

#Generate mesh
sim.mesh.genMesh(fnOut=fnMsh,debug=True)
	
	