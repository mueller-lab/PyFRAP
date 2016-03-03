import pyfrp_gmsh_IO_module 
import  pyfrp_term_module
import pyfrp_gmsh_geometry
import pyfrp_misc_module

meshFolder="meshfiles/"

files=pyfrp_misc_module.getSortedFileList(meshFolder,'geo')

for f in files:
	
	d,parmDic=pyfrp_gmsh_IO_module.readGeoFile(meshFolder+f)
	d.draw(ann=True)

raw_input()


#a,idx=d.getEdgeById(13)
#pyfrp_term_module.printAllObjAttr(a)

#pyfrp_term_module.printDict(parmDic)
