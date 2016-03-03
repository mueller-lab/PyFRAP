import sys
import meshio
		
fn=sys.argv[1]

#MeshIO

points, cells, point_data, cell_data, field_data = meshio.read(fn)
fnOut=fn.replace('.msh','.vtk')

meshio.write(fnOut,points,cells,point_data=point_data,cell_data=cell_data,field_data=field_data)