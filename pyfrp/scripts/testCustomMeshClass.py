###Script that creates custom mesh class and that runs Gmsh

#Import module
import pyfrp_gmsh_module

#Custom mesh class
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

m=mesh(True,fnGeo,fnMsh)

m.genMesh(fnOut=fnMsh,debug=True)

#Run
#pyfrp_gmsh_module.runGmsh(fnGeo,fnOut=fnMsh,debug=debug)

