#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFRAP.

#PyFRAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Mesh class for PyFRAP toolbox, including following classes:

#(1) mesh 

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================


#PyFRAP
from pyfrp.modules import pyfrp_gmsh_module
from pyfrp.modules import pyfrp_gmsh_IO_module
from pyfrp.modules import pyfrp_integration_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules.pyfrp_term_module import *

#FiPy
import fipy

#Numpy/Scipy
import numpy as np

#Misc
import os


#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================

class mesh(object):

	def __init__(self,simulation):
		
		#Naming/ID
		self.simulation=simulation
		self.mesh=None
		self.restoreDefaults()
	
	def setVolSizePx(self,v,remesh=True,fnOut=None):
		
		"""Sets volSize of mesh in px.
		
		.. note:: If ``fnOut=None``, then either ``fnMesh`` will be used, or,
		   if ``fnMesh`` is not set yet, will use ``geometry.fnGeo``.
		
		Args:
			v (float): New volSize.
			
		Keyword Args:
			remesh (bool): Generate mesh with new volSize.
			fnOut (str): Output filepath for meshfile.
		
		Returns:
			float: New volSize.
		
		"""
		
		
		self.volSizePx=v
		self.updateGeoFile()
		if remesh:
			self.genMesh(fnOut=fnOut)
		return self.volSizePx
	
	def getVolSizePx(self):
		
		"""Returns mesh volSize in px.
		
		Returns:
			float: VolSize.
		
		"""
		
		return self.volSizePx
		
	def getSimulation(self):
		
		"""Returns :py:class:`pyfrp.subclasses.pyfrp_simulation` that mesh belongs to.
		
		Returns:
			pyfrp.subclasses.pyfrp_simulation: Simulation object.
		"""
		
		return self.simulation
	
	def restoreDefaults(self):	
		
		"""Restores default parameters of mesh.
		
		Default parameters are:
		
			* ``mesh.geometry=mesh.simulation.embryo.geometry``
			* ``mesh.fromFile=True``
			* ``mesh.volSizePx=20``
			* ``mesh.fnMesh=""``
	
		"""
		
		if self.simulation.embryo.geometry==None:
			printWarning("Geometry of embryo not specified, mesh generation will not work!")
		
		self.geometry=self.simulation.embryo.geometry
		self.fromFile=True
		self.volSizePx=20
		self.fnMesh=""
		
	def genMesh(self,fnOut=None,debug=False):
		
		"""Main mesh generation function.
		
		.. note:: If ``fnOut=None``, will use ``geometry.fnGeo``.
		
		.. note:: If ``fromFile=True``, will generate from .geo file running
		   gmsh directly on the file. If not, will try to run hard coded FiPy
		   version for mesh generation via :py:func:`runFiPyMeshGenerator` .
		
		Keyword Args:
			fnOut (str): Output filepath for meshfile.
			debug (bool): Print debugging messages.
		
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
		
		"""
		
		if fnOut==None:
			fnOut=self.geometry.fnGeo.replace(".geo",".msh")
			
		if self.fromFile:
			self.fnMesh=fnOut
		
			pyfrp_gmsh_module.runGmsh(self.geometry.fnGeo,fnOut=fnOut,debug=debug,volSizeMax=self.volSizePx)
			
			self.importMeshFromFile(self.fnMesh)
		else:
			self.runFiPyMeshGenerator(self.geometry.typ)

		return self.mesh
	
	def runFiPyMeshGenerator(self,typ):
		
		"""Runs gmsh on the via FiPy internally defined meshes.
		
		Available meshes:
			
			* cylinder
			* zebrafishDomeStage
			* xenopusBall
			
		.. note:: Any refinement method will not work if mesh is created this way.
		
		Args:
			typ (str): Type of mesh to be created (see list above).
				
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
			
		"""
		
		if typ=="cylinder":
			self.mesh=pyfrp_gmsh_module.genFiPyCylinderMesh(self.volSizePx,self.geometry.radius,self.geometry.height,self.geometry.center)
		elif typ=="zebrafishDomeStage":
			self.mesh=pyfrp_gmsh_module.genFiPyDomeMesh(self.volSizePx,self.geometry.innerRadius,self.geometry.outerRadius,self.geometry.centerDist,self.geometry.center)
		elif typ=="xenopusBall":
			self.mesh=pyfrp_gmsh_module.genFiPyBallMesh(self.volSizePx,self.geometry.radius,self.geometry.center)
		else:
			printWarning("Geometry type " + typ + "unknown for runFiPyMeshGenerator. Check pyfrp_gmsh_module for available geometries.")
		
		return self.mesh
		
	def importMeshFromFile(self,fn):
		
		"""Imports mesh from a Gmsh .msh file.
		
		See also http://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html. 
		
		Args:
			fn (str): Filepath to meshfile.
				
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
			
		"""
		
		self.mesh=fipy.GmshImporter3D(fn)
		return self.mesh
		
	def setMesh(self,m):
		
		"""Sets ``mesh`` attribute to a new mesh.
		
		Args:
			m (fipy.GmshImporter3D): New mesh.
				
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
		
		"""
		
		self.mesh=m
		return self.mesh
	
	def setFromFile(self,v):
		
		"""Sets flag if mesh is supposed to be created from file (recommended)
		or from internally defined mesh creation method.
		
		Args:
			v (bool): New flag value.
				
		Returns:
			bool: New flag value.
		
		"""
		
		self.fromFile=v
		return self.fromFile
	
	def getMesh(self):
		
		"""Returns mesh that is used for simulation.
		
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
		
		"""

		return self.mesh
	
	def getFnMesh(self):
		
		"""Returns the filepath of meshfile.
		"""
		
		return self.fnMesh
	
	def setFnMesh(self,fn):
		
		
		"""Sets the filepath of meshfile.
		
		Imports the new mesh right away using :py:func:`importMeshFromFile`.
		
		"""
		
		
		self.fnMesh=fn
		self.importMeshFromFile(self.fnMesh)
		return self.fnMesh
	
	def updateGeoFile(self,debug=False):
		
		"""Updates geometry file by writing new ``volSizePx``
		into ``embryo.geometry.fnGeo``.
		
		Keyword Args:
			debug (bool): Print debugging messages.
		
		Returns:
			str: Path to ouput meshfile.
		
		"""
		
		if os.path.isfile(self.geometry.fnGeo):
			self.fnMesh=pyfrp_gmsh_module.updateVolSizeGeo(self.geometry.fnGeo,self.volSizePx,run=False,debug=debug)
		else:
			printWarning("Cannot update meshfile, since geometry.fnGeo does not exist or is invalid.")
		return self.fnMesh
	
	def refine(self,debug=False):
		
		"""Refines mesh by splitting.
		
		See also http://gmsh.info/doc/texinfo/gmsh.html .
		
		Keyword Args:
			debug (bool): Print debugging messages.
		
		"""
		
		pyfrp_gmsh_module.refineMsh(self.fnMesh,debug=debug)
		self.importMeshFromFile(self.fnMesh)
		return self.fnMesh
	
	def forceMinMeshDensityInROI(self,ROI,density,stepPercentage=0.1,debug=False,findIdxs=True,method='refine',maxCells=100000):
		
		"""Forces global mensh density such that a certain density is reached in a
		given ROI.
		
		Tries to achive a mesh density ``density`` in ``ROI`` by globally refining mesh either 
		through decreasing ``volSizePx`` by ``stepPercentage`` percent (``method=volSize``), or 
		by using Gmsh's ``-refine`` option (``method=refine``). If maximum number of cells is 
		exceeded, will use the last mesh that did not exceed ``maxCells``.
		
		Args:
			ROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI object.
			density (float): Desired density.
			
		Keyword Args:
			stepPercentage (float): If method is ``volSize``, percentage of ``volSize`` decrease.
			method (str): Refinement method (``refine``/``volSize``).
			maxCells (int): Total maximum number of mesh cells allowed.
			findIdxs (bool): Find ROI indices after refinement.
			debug (bool): Print debugging messages.
		
		Returns:
			float: New ``volSizePx``
		
		"""
		
		#Set counter to 0
		j=0
		
		#Loop until mesh density has been found
		while ROI.getMeshDensity()<density:
			
			#Backup old stuff
			volSizeBackup=self.getVolSizePx()
			meshBackup=self.getMesh()
			
			#Refine
			if method=='volSize':
				self.setVolSizePx((1-stepPercentage)*self.getVolSizePx())
			elif method=='refine':
				self.refine()
			else:
				printError("Unknown method: ", method)
				return
			
			#Recompute Idxs
			ROI.computeMeshIdx(self.mesh)
			
			#Debugging output
			if debug:
				if method=='volSize':
					print "Tried volSizePx ", self.getVolSizePx()
				elif method=='refine':
					print "Refinement step, ", j
				
				print "Current density ", ROI.getMeshDensity(), " desired density ", density
				print "Current number of nodes in ROI ", len(ROI.meshIdx)
				print " Mesh now has ", np.shape(self.mesh.x)[0], "cells." 
			
			#Check if we reached maxCells
			if np.shape(self.mesh.x)[0]>maxCells:
				if debug:
					printWarning("Reached maxCells, will keep previous iteration.")
				
				self.mesh=meshBackup
				self.setVolSizePx(volSizeBackup,remesh=False)
				
			#Increment counter
			j=j+1
			
		if debug:
			print "volSizePx = ", self.getVolSizePx(), "is sufficient."
		
		#Recompute idxs for all ROIs
		if findIdxs:
			self.simulation.embryo.computeROIIdxs()
					
		return self.getVolSizePx()	
	
	def getNNodes(self):
		
		"""Returns number of nodes in mesh. 
		
		If no mesh has been generated yet, will return 0.
		
		Returns:
			int: Number of nodes.
		
		"""
		
		if self.mesh==None:
			return 0
		else:
			return len(self.mesh.getCellCenters()[0])
		
	def writeVTKFile(self,fn="",sub=False):
		
		"""Writes mesh into vtk file.
		
		Uses *meshIO* (https://github.com/nschloe/meshio), to convert the mesh saved
		in ``fnMesh`` to a .vtk file. 
		
		If ``sub==True``, will start a seperate subprocess and submit 
		``pyfrp_meshIO_script.py`` to it. This can be sometimes useful, since PyFRAP
		sometimes tends to crash otherwise.
		
		If no output path is given via ``fn``, will use same path as ``fnMesh``.
		
		.. note:: *meshIO* only gets imported inside this function, making PyFRAP running
		   even without the package installed. However, this feature will only run with
		   *meshIO*.
		
		Keyword Args:
			fn (str): Optional output path.
			sub (bool): Subprocess flag.
		
		Returns:
			str: Used output path.
		
		"""
		
		
		if not os.path.isfile(self.fnMesh):
			printWarning("Filepath to meshfile has not been specified yet. Cannot write VTK file.")
		
		if fn=="":
			fn=self.fnMesh.replace('.msh','.vtk')
		
		if sub:
			
			cmd = "python pyfrp_meshIO_script.py "+ self.fnMesh
			
			import shlex
			import subprocess
			args = shlex.split(cmd)
			
			p = subprocess.Popen(args)
			p.wait()
		else:
			#MeshIO
			import meshio
			
			points, cells, point_data, cell_data, field_data = meshio.read(self.fnMesh)
			meshio.write(fn,points,cells,point_data=point_data,cell_data=cell_data,field_data=field_data)
			
		return fn
	
	def importVTKFile(self,fnVTK="",sub=False):
		
		"""Imports a .vtk file into a vtk renderer.
		
		If ``fnVTK`` is not given, will generate .vtk file from meshfile stored in ``fnMesh``
		using :py:func:`writeVTKFile`.
		
		If ``sub==True``, will start a seperate subprocess and submit 
		``pyfrp_meshIO_script.py`` to it. This can be sometimes useful, since PyFRAP
		sometimes tends to crash otherwise.
		
		.. note:: This function imports *vtk*. *vtk* is only necessary in a few functions,
		   hence only imported when needed. This should make PyFRAP more portable.
		
	
		Keyword Args:
			fnVTK (str): Path to input vtk file.
			sub (bool): Subprocess flag.
			
		Returns:
			vtk.vtkRenderer: Renderer object.
		
		"""
		
		#vtk
		import vtk
		
		if not os.path.isfile(self.fnMesh):
			printWarning("Filepath to meshfile has not been specified yet. Cannot plot.")
		
		if fnVTK=="":
			fnVTK=self.writeVTKFile(sub=sub)
		
		# Read the source file.
		reader = vtk.vtkUnstructuredGridReader()
		reader.SetFileName(fnVTK)
		reader.Update() 
		output = reader.GetOutput()
		
		#Extract Edges 
		edges=vtk.vtkExtractEdges()
		edges.SetInput(reader.GetOutput()) 
		
		#Make edges into tubes
		tubes = vtk.vtkTubeFilter()
		tubes.SetInput(edges.GetOutput())
		tubes.SetRadius(0.5)
		tubes.SetNumberOfSides(3)

		#Genereate wireframe mapper
		wireFrameMapper=vtk.vtkPolyDataMapper()
		wireFrameMapper.SetInput(tubes.GetOutput())
		wireFrameMapper.SetScalarVisibility(0)
		
		#Make Actor
		wireFrameActor=vtk.vtkActor()
		wireFrameActor.SetMapper(wireFrameMapper)
		wireFrameActor.GetProperty().SetColor(0,0,0)
		wireFrameActor.SetPickable(0)

		#Create the Renderer
		renderer = vtk.vtkRenderer()
		renderer.AddActor(wireFrameActor)
		renderer.SetBackground(1, 1, 1) # Set background to white

		return renderer
	
	def plotMesh(self,fnVTK=""):
		
		"""Plots the mesh using VTK.
		
		If ``fnVTK`` is not given, will generate .vtk file from meshfile stored in ``fnMesh``
		using :py:func:`writeVTKFile`.
		
		.. note:: This function imports *vtk*. *vtk* is only necessary in a few functions,
		   hence only imported when needed. This should make PyFRAP more portable.
		
		Keyword Args:
			fnVTK (str): Path to input vtk file.
			
		Returns:
			vtk.vtkRenderWindow: RenderWindow object.
		
		"""
		
		
		#vtk
		import vtk
		
		#import vtk file
		renderer=self.importVTKFile(fnVTK=fnVTK)
		
		#Window
		renderWindow = vtk.vtkRenderWindow()
		renderWindow.AddRenderer(renderer)
		
		#Interactor
		interactor = vtk.vtkRenderWindowInteractor()
		interactor.SetRenderWindow(renderWindow)
		
		#Start
		renderWindow.GetInteractor().Initialize()
		renderWindow.GetInteractor().Start()
		
		return renderWindow
		
	def printStats(self,tetLenghts=False):
		
		"""Prints out statistics of mesh.
		
		Also calculates all tetraheder lengths if ``tetLenghts`` is selected. 
		This might take some time depending on mesh size.
		
		Keyword Args:
			tetLenghts (bool): Also calculate and print out tetrahedra sidelengths.
		
		"""
	
			
		print "-------------------------------------------"
		print "Mesh Statistics:"
		print "-------------------------------------------"
		print "Mesh has ", np.shape(self.mesh.x)[0] , " cells"
		print "Mesh has ", self.mesh._numberOfVertices, " vertices"
		print "Mesh has ", self.mesh.numberOfFaces, " faces"
		print "min x=", min(self.mesh.x), "max x=", max(self.mesh.x)
		print "min y=", min(self.mesh.y), "max y=", max(self.mesh.y)
		print "min z=", min(self.mesh.z), "max z=", max(self.mesh.z)
		print "Maximum cell volume= ", max(self.mesh.getCellVolumes())
		print "Minimum cell volume= ", min(self.mesh.getCellVolumes())
			
		print "Maximum cell volume is", max(self.mesh.getCellVolumes()), "in cell number=", np.argmax(self.mesh.getCellVolumes())
		print "Minimum cell volume is", min(self.mesh.getCellVolumes()), "in cell number=", np.argmin(self.mesh.getCellVolumes())
		print "Average cell volume is", np.mean(self.mesh.getCellVolumes())

		if tetLenghts:
			slsVec=self.calcAllTetSidelenghts()
			print "Average sidelength of tetrahedron in self.mesh:", np.mean(slsVec)
			print "Maximum sidelength of tetrahedron in self.mesh:", max(slsVec)
			print "Minimum sidelength of tetrahedron in self.mesh:", min(slsVec)	
		
		print
		
		return
	
	def calcAllTetSidelenghts(self):
		
		"""Calculates sidelengths of all tetrahedra.
		
		See also :py:func:`pyfrp.modules.pyfrp_integration_module.calcTetSidelengths`.
		
		Returns:
			list: List of all sidelengths.
		"""
		
		#Calculating sidelengths of tetrahedra
		slsVec=[]
		for i in range(np.shape(self.mesh._getOrderedCellVertexIDs())[1]):
			currVert=self.mesh._getOrderedCellVertexIDs()[:,i]
			
			pt1=self.mesh.vertexCoords[:,currVert[0]]
			pt2=self.mesh.vertexCoords[:,currVert[1]]
			pt3=self.mesh.vertexCoords[:,currVert[2]]
			pt4=self.mesh.vertexCoords[:,currVert[3]]
			
			sl1,sl2,sl3=pyfrp_integration_module.calcTetSidelengths(pt1,pt2,pt3,pt4)
			
			slsVec.append(sl1)
			slsVec.append(sl2)
			slsVec.append(sl3)
			
		return slsVec
	
	def plotDensity(self,axes=None,hist=True,bins=100,color='b'):
		
		"""Plots the mesh density in x/y/z-direction.
		
		``hist=True`` is recommended, since otherwise plots generally appear fairly 
		noisy. 
		
		.. note:: If no ``axes`` are given or they do not have the necessary size, 
		   will create new ones.
		
		.. image:: ../imgs/pyfrp_mesh/density_plot.png
		
		Keyword Args:
			axes (list): List of ``matplotlib.axes``.
			hist (bool): Summarize densities in bins.
			bins (int): Number of bins used for hist.
			color (str): Color of plot.
			
		Returns:
			list: List of ``matplotlib.axes``.
			
		"""
			
		volSortedByX,xSorted=pyfrp_misc_module.sortListsWithKey(self.mesh.getCellVolumes(),self.mesh.x)
		volSortedByY,ySorted=pyfrp_misc_module.sortListsWithKey(self.mesh.getCellVolumes(),self.mesh.y)
		volSortedByZ,zSorted=pyfrp_misc_module.sortListsWithKey(self.mesh.getCellVolumes(),self.mesh.z)
		
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["Density(x)","Density(y)","Density(z)"])
		else:
			if len(axes)<3:
				printWarning("axes do not have right have, will create new ones.")
				fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["Density(x)","Density(y)","Density(z)"])
		
		if hist:
			xSorted,volSortedByX=pyfrp_misc_module.simpleHist(xSorted,volSortedByX,bins)
			ySorted,volSortedByY=pyfrp_misc_module.simpleHist(ySorted,volSortedByY,bins)
			zSorted,volSortedByZ=pyfrp_misc_module.simpleHist(zSorted,volSortedByZ,bins)
			
		axes[0].plot(xSorted,volSortedByX,color=color)
		axes[1].plot(ySorted,volSortedByY,color=color)
		axes[2].plot(zSorted,volSortedByZ,color=color)
		
		for ax in axes:
			pyfrp_plot_module.redraw(ax)
		
		return axes
			
	def addBoxField(self,volSizeIn,rangeX,rangeY,rangeZ,newFile=True,fnAppendix="_box",comment="newField",run=False):
		
		"""Adds box field to mesh.
		
		Box fields allow to refine certain areas of the mesh, see also 
		http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes .
		
		.. note:: Will keep ``volSizePx`` as volSize outside of the box.
		
		Args:
			volSizeIn (float): volSize in px inside the box.
			rangeX (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
			rangeY (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
			rangeZ (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
			newFile (bool): Write new mesh into a new .geo file.
			fnAppendix (str): Append this to new file name.
			comment (str): Comment in .geo file before definition of box field.
			run (bool): Run Gmsh on new .geo file afterwards.
			
		Returns:
			str: Path to new .geo file.
			
		"""
		
		if newFile:
			if "_box" not in self.geometry.fnGeo:	
				fnOut=os.path.dirname(self.geometry.fnGeo)+"/field/custom/"+os.path.basename(self.geometry.fnGeo).replace(".geo",fnAppendix+"_"+self.geometry.embryo.name+".geo")
			else:
				fnOut=fnOut=self.geometry.fnGeo
		else:
			fnOut=self.geometry.fnGeo
		
		pyfrp_gmsh_IO_module.addBoxField(self.geometry.fnGeo,volSizeIn,self.volSizePx,rangeX,rangeY,rangeZ,comment=comment,fnOut=fnOut)
		self.geometry.setFnGeo(fnOut)
		
		
		if run:
			self.genMesh()
		
		return fnOut
	
			
		
		
		