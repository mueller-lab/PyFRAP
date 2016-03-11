###Script to debug mesh class on cluster
#Idea: Generate everything is normal, but then sub sim.mesh with minimal mesh class

#Import module
import time
import pyfrp_embryo
import pyfrp_gmsh_module

#Define minimal mesh class
class mesh(object):
	
	def __init__(self,fromFile,fnGeo,fnMsh):
		
		self.fromFile=fromFile
		self.fnGeo=fnGeo
		self.fnMsh=fnMsh

	def genMesh(self,fnOut=None,debug=False):
		if fnOut==None:
			fnOut=self.fnGeo.replace(".geo",".msh")
		
		if self.fromFile:
			self.fnMsh=fnOut
			pyfrp_gmsh_module.runGmsh(self.fnGeo,fnOut=fnOut,debug=debug)

#Define variables
fnGeo='meshfiles/cylinder.geo'
fnMsh='meshfiles/cylinder.msh'
debug=True

#Geometrical parameters
#center=[256,256]
#cylinderHeight=90.3332804037
#imagingRadius=294.103867936	
	
##Create embryo	
#emb=pyfrp_embryo.embryo("TestedOnCluster")	

##Create geometry
#emb.setGeometry2Cylinder(center,imagingRadius,cylinderHeight)
	
##Create Simulation
#sim=emb.newSimulation()	

#Now substitute with minimal mesh
m=mesh(True,fnGeo,fnMsh)
m.genMesh(fnOut=fnMsh,debug=True)
#sim.mesh.setMesh(m)

#Run Gmsh

	
	
#m=mesh(True,fnGeo,fnMsh)

#m.genMesh(fnOut=fnMsh,debug=True)