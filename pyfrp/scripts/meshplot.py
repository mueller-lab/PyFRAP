#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

import pyfrp_IO_module
import sys
from pyfrp_term_module import *
import numpy as np
from vtk import *

import meshio

#===========================================================================================================================================================================
#Run script
#===========================================================================================================================================================================

#Pass filename
fn=sys.argv[1]

#Load Embryo
mol=pyfrp_IO_module.loadFromPickle(fn)

emb=mol.embryos[-1]

print emb

m=emb.simulation.mesh

m.plotMesh()

raw_input()

#ordvert=m._getOrderedCellVertexIDs()

#print m

#print ordvert

#dir(m.x)

#raw_input()

#points, cells, point_data, cell_data, field_data = meshio.read('meshfiles/dome.msh')

#meshio.write('meshfiles/dome.vtk',points,cells,point_data=point_data,cell_data=cell_data,field_data=field_data)

## The source file
#file_name = "meshfiles/dome.vtk"
 
## Read the source file.
#reader = vtkUnstructuredGridReader()
#reader.SetFileName(file_name)
#reader.Update() # Needed because of GetScalarRange
#output = reader.GetOutput()
#scalar_range = output.GetScalarRange()


##Extract Edges 
#edges=vtkExtractEdges()
#edges.SetInput(reader.GetOutput()) 
 
#tubes = vtkTubeFilter()
#tubes.SetInput(edges.GetOutput())
#tubes.SetRadius(0.5)
#tubes.SetNumberOfSides(3)

#wireFrameMapper=vtkPolyDataMapper()
#wireFrameMapper.SetInput(tubes.GetOutput())
#wireFrameMapper.SetScalarVisibility(0)

#wireFrameActor=vtkActor()
#wireFrameActor.SetMapper(wireFrameMapper)
#wireFrameActor.GetProperty().SetColor(0,0,0)
#wireFrameActor.SetPickable(0)

 
## Create the mapper that corresponds the objects of the vtk file
## into graphics elements
##mapper = vtkDataSetMapper()
##mapper.SetInput(output)
##mapper.SetScalarRange(scalar_range)
##mapper.ScalarVisibilityOff() 
 
## Create the Actor
##actor = vtkActor()
##actor.SetMapper(mapper)
##actor.GetProperty().SetRepresentationToWireframe()
##actor.GetProperty().SetColor(0,1,0)


##polydatamapper
 
## Create the Renderer
#renderer = vtkRenderer()
#renderer.AddActor(wireFrameActor)
#renderer.SetBackground(1, 1, 1) # Set background to white

 
## Create the RendererWindow
#renderer_window = vtkRenderWindow()
#renderer_window.AddRenderer(renderer)
 
## Create the RendererWindowInteractor and display the vtk_file
#interactor = vtkRenderWindowInteractor()
#interactor.SetRenderWindow(renderer_window)
#interactor.Initialize()
#interactor.Start()





		#// mesh is my unstructuredgrid

		#surfFilter=new vtkDataSetSurfaceFilter();
		#surfFilter.SetInput(mesh);

		#vtkExtractEdges edges = new vtkExtractEdges();
		#edges.SetInput(surfFilter.GetOutput());
		
		#vtkTubeFilter tubes = new vtkTubeFilter();
		#tubes.SetInput(edges.GetOutput());
		#tubes.SetRadius(tubeRadius);
		#tubes.SetNumberOfSides(10);
		

		#wireFrameMapper=new vtkPolyDataMapper();
		#wireFrameMapper.SetInput(tubes.GetOutput());
		#wireFrameMapper.SetScalarVisibility(0);	
		
		#wireFrameActor=new vtkActor();
		#wireFrameActor.SetMapper(wireFrameMapper);
		#wireFrameActor.GetProperty().SetColor(0,0,0);
		#wireFrameActor.SetPickable(0);   

		#canvas.GetRenderer().AddActor(wireFrameActor);