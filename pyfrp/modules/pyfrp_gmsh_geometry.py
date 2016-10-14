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

"""PyFRAP module for creating/extracting gmsh geometries for PyFRAP toolbox. Module mainly has the following classes:

	* A ``domain`` class, acting as a canvas.
	* A ``vertex`` class, substituting gmsh's *Point*.
	* A ``edge`` class, parenting all different kind of edges.
	* A ``line`` class, substituting gmsh's *Line*.
	* A ``arc`` class, substituting gmsh's *Circle*.
	
This module together with pyfrp.pyfrp_gmsh_IO_module and pyfrp.pyfrp_gmsh_module works partially as a python gmsh wrapper, however is incomplete.
If you want to know more about gmsh, go to http://gmsh.info/doc/texinfo/gmsh.html .
"""
#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#String
import string

#PyFRAP modules
import pyfrp_plot_module
from pyfrp_term_module import *
import pyfrp_misc_module
import pyfrp_gmsh_IO_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def getAngle(vec1,vec2):
	
		"""Returns angle between two vectors in radians.
		
		Args:
			vec1 (numpy.ndarray): Vector 1.
			vec2 (numpy.ndarray): Vector 2.
		
		Returns: 
			float: Angle.
		
		"""
		
		a=np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
		
		if a<0:
			return getAngle(vec2,vec1)
		return a

def flipCoordinate(x,destAxis,origAxis='x',debug=False):
	
		"""Transforms coodinate from one axis to another by
		rolling the coordinates, e.g. clockwise turning the 
		point.
		
		``destAxis`` and ``origAxis`` are given as one of 
		``x,y,z``. 
		
		Args:
			x (numpy.ndarray): Coordinate to turn.
			destAxis (str): Destination axis.
			
		Keyword Args:	
			origAxis (str): Original axis.
			debug (bool): Print debugging output.
		
		Returns:
			numpy.ndarray: Transformed coordinate.
		
		"""
		
			
		# Calculate differences between axis
		axisDiff=abs(string.lowercase.index(destAxis)-string.lowercase.index(origAxis))
		
		# Roll
		xnew=np.roll(x,axisDiff)
		
		# Print debugging messages
		if debug:
			print "Transforming coordinate " , x, " from axis ", origAxis, " to axis ", destAxis , "."
			print "axisDiff = ", axisDiff
			print "xnew = ", xnew
		
		return xnew 
		
#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================
	
class domain:
	
	"""Domain class storing embryo geometry entities.
	
	Args:
		edges (list): List of edges.
		vertices (list): List of vertices.
		arcs (list): List of arcs.
		lines (list): List of lines.
		lineLoops (list): List of lineLoops.
		surfaceLoops (list): List of surfaceLoops.
		ruledSurfaces (list): List of ruledSurfaces.
		volumes (list): List of volumes.
		annXOffset (float): Offset of annotations in x-direction.
		annYOffset (float): Offset of annotations in y-direction.
		annZOffset (float): Offset of annotations in z-direction.
			
	"""
	
	def __init__(self):
		
		self.edges=[]
		self.vertices=[]
		self.arcs=[]
		self.lines=[]
		self.lineLoops=[]
		self.ruledSurfaces=[]
		self.surfaceLoops=[]
		self.volumes=[]
		
		self.annXOffset=3.
		self.annYOffset=3.
		self.annZOffset=3.
		
		
	def addVertex(self,x,Id=None,volSize=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` instance
		at point ``x`` and appends it to ``vertices`` list.
		
		.. note:: ``volSize`` does not have any effect on the geometry itself but is simply 
		   stored in the vertex object for further usage.
		
		Args:
			x (numpy.ndarray): Coordinate of vertex.
			
		Keyword Args:
			Id (int): ID of vertex.
			volSize (float): Element size at vertex.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.vertex: New vertex instance.
		
		"""
		
		newId=self.getNewId(self.vertices,Id)
		
		v=vertex(self,x,newId,volSize=volSize)
		self.vertices.append(v)	
		
		return v
	
	def addLine(self,v1,v2,Id=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.line` instance
		at point ``x`` and appends it to ``edges`` and ``lines`` list.
		
		Args:
			v1 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
			v2 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
			
		Keyword Args:
			Id (int): ID of line.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.line: New line instance.
		
		"""
		
		newId=self.getNewId(self.lines,Id)
		
		e=line(self,v1,v2,newId)
		self.lines.append(e)
		self.edges.append(e)
		
		return e
	
	def addArc(self,vstart,vcenter,vend,Id=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.arc` instance
		at point ``x`` and appends it to ``edges`` and ``arcs`` list.
		
		Args:
			vstart (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
			vcenter (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Center vertex.
			vend (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
			
		Keyword Args:
			Id (int): ID of arc.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.arc: New line instance.
		
		"""
		
		
		newId=self.getNewId(self.arcs,Id)
			
		a=arc(self,vstart,vcenter,vend,newId)
		self.arcs.append(a)
		self.edges.append(a)
		
		return a
	
	def addCircleByParameters(self,center,radius,z,volSize,plane="z"):
		
		"""Adds circle to domain by given center and radius.
		
		Will create 5 new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` objects 
		``[vcenter,v1,v2,v3,v4]`` and four new `pyfrp.modules.pyfrp_gmsh_geometry.arc` objects
		[a1,a2,a3,a4] and builds circle.
		
		Circle  will be at ``z=z`` and vertices will have mesh size ``volSize``.
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addCircleByParameters.png
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.flipCoordinate`.
		
		Args:
			center (numpy.ndarray): Center of circle.
			radius (float): Radius of the circle.
			z (float): Height at which circle is placed.
			volSize (float): Height at which circle is placed.
		
		Keyword Args:
			plane (str): Plane in which circle is placed.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.arc: New line instance.
		
		"""
		
		# Define coordinates
		xcenter=flipCoordinate([center[0],center[1],z],plane,origAxis="z")
		x1=flipCoordinate([center[0]+radius,center[1],z],plane,origAxis="z")
		x2=flipCoordinate([center[0],center[1]+radius,z],plane,origAxis="z")
		x3=flipCoordinate([center[0]-radius,center[1],z],plane,origAxis="z")
		x4=flipCoordinate([center[0],center[1]-radius,z],plane,origAxis="z")
		
		# Add vertices
		vcenter=self.addVertex(xcenter,volSize=volSize)
		v1=self.addVertex(x1,volSize=volSize)
		v2=self.addVertex(x2,volSize=volSize)
		v3=self.addVertex(x3,volSize=volSize)
		v4=self.addVertex(x4,volSize=volSize)
		
		# Add Arcs
		a1=self.addArc(v1,vcenter,v2)
		a2=self.addArc(v2,vcenter,v3)
		a3=self.addArc(v3,vcenter,v4)
		a4=self.addArc(v4,vcenter,v1)
		
		return [vcenter,v1,v2,v3,v4],[a1,a2,a3,a4]
		
		
	#def addCylinderByParamters(self,center,radius,z,height,volSize,plane="z"):
		
		###NOTE: CONTINUE here
		
	def checkIdExists(self,Id,objList):
		
		"""Checks if any object in ``objList`` already has ID ``Id``.
		
		Args:
			Id (int): ID to be checked.
			objList (list): List of objects, for example ``edges``.
			
		Returns:
			bool: True if any object has ID ``Id``.
		
		"""
		
		IdList=pyfrp_misc_module.objAttrToList(objList,'Id')
		if Id in IdList:
			printWarning("Object with Id " + str(Id) + " already exists.")
			return True
		return False
	
	def getNewId(self,objList,Id=None):
		
		"""Returns free ID for object type.
		
		Args:
			objList (list): List of objects, for example ``edges``.
			
		Keyword Args:
			Id (int): ID to be checked.
			
		Returns:
			int: New free ID.
		
		"""
		
		if Id==None:
			newId=self.incrementID(objList)
		else:
			if self.checkIdExists(Id,objList):
				newId=self.incrementID(objList)
			else:
				newId=Id
		
		return newId
		
	def incrementID(self,objList):
		
		"""Returns ID that is by one larger for a specific 
		object type.
		
		Args:
			objList (list): List of objects, for example ``edges``.
				
		Returns:
			int: Incremented ID.
		
		"""
		
		if len(objList)==0:
			newId=0
		else:
			IdList=pyfrp_misc_module.objAttrToList(objList,'Id')
			newId=max(IdList)+1		
		return newId
		
	def getEdgeById(self,ID):
		
		"""Returns edge with ID ``ID``.
		
		Returns ``(False,False)`` if edge cannot be found.
		
		Args:
			ID (int): ID of edge.
				
		Returns:
			tuple: Tuple containing:
				
				* e (pyfrp.modules.pyfrp_gmsh_geometry.edge): Edge.
				* i (int): Position in ``edges`` list.
		
		"""
		
		for i,e in enumerate(self.edges):
			if e.Id==ID:
				return e,i
		return False,False
	
	def getLineLoopById(self,ID):
		
		"""Returns lineLoop with ID ``ID``.
		
		Returns ``(False,False)`` if lineLoop cannot be found.
		
		Args:
			ID (int): ID of lineLoop.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): lineLoop.
				* i (int): Position in ``lineLoops`` list.
		
		"""
		
		for i,l in enumerate(self.lineLoops):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getRuledSurfaceById(self,ID):
		
		"""Returns ruledSurface with ID ``ID``.
		
		Returns ``(False,False)`` if ruledSurface cannot be found.
		
		Args:
			ID (int): ID of ruledSurface.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): ruledSurface.
				* i (int): Position in ``ruledSurfaces`` list.
		
		"""
		
		for i,l in enumerate(self.ruledSurfaces):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getSurfaceLoopById(self,ID):
		
		"""Returns surfaceLoop with ID ``ID``.
		
		Returns ``(False,False)`` if surfaceLoop cannot be found.
		
		Args:
			ID (int): ID of surfaceLoop.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): surfaceLoop.
				* i (int): Position in ``surfaceLoops`` list.
		
		"""
		
		for i,l in enumerate(self.surfaceLoops):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getVolumeById(self,ID):
		
		"""Returns volume with ID ``ID``.
		
		Returns ``(False,False)`` if volume cannot be found.
		
		Args:
			ID (int): ID of volume.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.volume): volume.
				* i (int): Position in ``volumes`` list.
		
		"""
		
		for i,l in enumerate(self.volumes):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getVertexById(self,ID):
		
		"""Returns vertex with ID ``ID``.
		
		Returns ``(False,False)`` if vertex cannot be found.
		
		Args:
			ID (int): ID of vertex.
				
		Returns:
			tuple: Tuple containing:
				
				* v (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Vertex.
				* i (int): Position in ``edges`` list.
		
		"""
		
		for i,v in enumerate(self.vertices):
			if v.Id==ID:
				return v,i
		return False,False
		
	def draw(self,ax=None,color=None,ann=None):
		
		"""Draws complete domain.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of domain.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if color==None:
			color='k'
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],proj=['3d'])
			ax=axes[0]
		
		for v in self.vertices:
			v.draw(ax=ax,color=color,ann=ann)
		for e in self.edges:
			e.draw(ax=ax,color=color,ann=ann)	
		for a in self.arcs:
			a.draw(ax=ax,color=color,ann=ann)	
		
		return ax
		
	def getExtend(self):
		
		"""Returns extend of domain in all 3 dimensions.
		
		Returns: 
			tuple: Tuple containing:
				
				* minx (float): Minimal x-coordinate.
				* maxx (float): Maximal x-coordinate.
				* miny (float): Minimal y-coordinate.
				* maxy (float): Maximal y-coordinate.
				* minz (float): Minimal z-coordinate.
				* maxz (float): Maximal z-coordinate.
	
		"""
			
		x=[]
		y=[]
		z=[]
		for v in self.vertices:
			x.append(v.x[0])
			y.append(v.x[1])
			z.append(v.x[2])
		return min(x), max(x), min(y),max(y), min(z),max(z)
	
	def verticesCoordsToList(self):
		
		"""Returns list of coordinates from all vertrices.
		
		Returns:
			list: List of (x,y,z) coordinates.
		
		"""
		
		l=[]
		for v in self.vertices:
			l.append(v.x)
		return l
	
	def addLineLoop(self,Id,edgeIDs=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.lineLoop` instance
		with given edgeIDs. 
			
		Keyword Args:
			edgeIDs (list): List of edge IDs included in line loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.lineLoop: New lineLoop instance.
		
		"""
		
		newId=self.getNewId(self.lineLoops,Id)
		
		l=lineLoop(self,edgeIDs,newId)
		self.lineLoops.append(l)
		
		return l
	
	def addSurfaceLoop(self,Id,surfaceIDs=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop` instance
		with given surfaceIDs. 
			
		Keyword Args:
			surfaceIDs (list): List of surface IDs included in surface loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop: New surfaceLoop instance.
		
		"""
		
		newId=self.getNewId(self.surfaceLoops,Id)
		
		l=surfaceLoop(self,surfaceIDs,newId)
		self.surfaceLoops.append(l)
		
		return l
	
	def addRuledSurface(self,Id,lineLoopID=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface` instance
		with given lineLoop. 
			
		Keyword Args:
			lineLoopID (ID): ID of line loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface: New ruledSurface instance.
		
		"""
		
		newId=self.getNewId(self.surfaceLoops,Id)
		
		l=ruledSurface(self,lineLoopID,newId)
		self.ruledSurfaces.append(l)
		
		return l
	
	def addVolume(self,Id,surfaceLoopID=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.volume` instance
		with given surfaceLoop. 
			
		Keyword Args:
			surfaceLoopID (ID): ID of surface loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.volume: New volume instance.
		
		"""
		
		newId=self.getNewId(self.volumes,Id)
		
		l=volume(self,surfaceLoopID,newId)
		self.volumes.append(l)
		
		return l
	
	
	def writeToFile(self,fn):
		
		"""Writes domain to file.
		
		Args:
			fn (str): File path to write to.
			
		"""
		
		with open(fn,'wb') as f:
			
			self.writeElements("vertices",f)
			self.writeElements("lines",f)
			self.writeElements("arcs",f)
			self.writeElements("lineLoops",f)
			self.writeElements("ruledSurfaces",f)
			self.writeElements("surfaceLoops",f)
			self.writeElements("volumes",f)

	def writeElements(self,element,f):
			
		"""Writes all entities of a specific element type to file.
		
		Possible elements are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* ruledSurfaces
			* surfaceLoops
			* volumes
		
		Args:
			element (str): Element type to write.
			f (file): File to write to.
			
		"""
	
		pyfrp_gmsh_IO_module.writeComment(f,element)
		for v in getattr(self,element):
			f=v.writeToFile(f)
		f.write("\n")
	
	def incrementAllIDs(self,offset):
		
		"""Adds offset to all entity IDs.
		
		Args:
			offset (int): Offset to be added.
			
		"""
		
		self.incrementIDs(offset,"vertices")
		self.incrementIDs(offset,"lines")
		self.incrementIDs(offset,"arcs")
		self.incrementIDs(offset,"lineLoops")
		self.incrementIDs(offset,"ruledSurfaces")
		self.incrementIDs(offset,"surfaceLoops")
		self.incrementIDs(offset,"volumes")
		
	def incrementIDs(self,offset,element):
		
		"""Adds offset to all entity IDs.
		
		Possible elements are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* ruledSurfaces
			* surfaceLoops
			* volumes
		
		Args:
			offset (int): Offset to be added.
			element (str): Element type to increment.
		
		"""
		
		for e in getattr(self,element):		
			e.Id=e.Id+offset
			
class vertex:
	
	"""Vertex class storing information from gmsh .geo Points.
	
	.. note:: ``volSize`` does not have any effect on the geometry itself but is simply 
		stored in the vertex object for further usage.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain vertex belongs to.
		x (numpy.ndarray): Coordinate of vertex.
		Id (int): ID of vertex.
			
	Keyword Args:
		volSize (float): Element size at vertex.
		
	"""
	
	def __init__(self,domain,x,Id,volSize=None):
		self.domain=domain
		
		self.x=np.array(x)
		self.Id=Id
		self.volSize=volSize
			
	def draw(self,ax=None,color=None,ann=None):
		
		"""Draws vertrex into axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of domain.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],proj=['3d'])
			ax=axes[0]
			
		ax.scatter(self.x[0],self.x[1],self.x[2],c=color)
		if ann:
			ax.text(self.x[0]+self.domain.annXOffset, self.x[1]+self.domain.annYOffset, self.x[2]+self.domain.annZOffset, "p"+str(self.Id), None)
			
		pyfrp_plot_module.redraw(ax)
		
		return ax
		
	def setX(self,x):
		
		"""Sets coordinate if vertex to ``x``.
		
		Returns:
			numpy.ndarray: New vertex coordinate.
		"""
		
		self.x=x
		return self.x
	
	def writeToFile(self,f):
		
		"""Writes vertex to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Point("+str(self.Id)+")= {" + str(self.x[0]) + ","+ str(self.x[1])+ "," + str(self.x[2]) + ',' + str(self.volSize) + "};\n" )
		
		return f
	
	
class edge:
	
	"""Edge class storing information from gmsh .geo circles and lines.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain vertex belongs to.
		Id (int): ID of edge.
		typ (int): Type of edge (1=arc/0=line).
			
	"""
	
	def __init__(self,domain,Id,typ):
		self.domain=domain
		self.Id=Id
		self.typ=typ
		
	def getDomain(self):
		
		"""Returns domain edge belongs to."""
		
		return self.domain
	
	def getId(self):
		
		"""Returns edges ID."""
		
		return self.Id
	
	def getTyp(self):
		
		"""Returns Type of edge."""
		
		return self.typ
	
	def decodeTyp(self):
		
		"""Decodes type of edge into string."""
		
		if typ==1:
			return "arc"
		elif typ==0:
			return "line"
	
class line(edge):
	
		
	"""Line class storing information from gmsh .geo lines.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain line belongs to.
		v1 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
		v2 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
		Id (int): ID of line.
		
	"""
	
	
	def __init__(self,domain,v1,v2,Id):
		
		edge.__init__(self,domain,Id,0)

		self.v1=v1
		self.v2=v2
	
	def getMiddle(self):
		
		r"""Returns midpoint of line.
		
		.. math:: m = \frac{x(v_1) + x(v_2)}{2}
		     
		Returns:
			numpy.ndarray: Midpoint.
			
		"""
		
		return (self.v1.x+self.v2.x)/2.
	
	def draw(self,ax=None,color=None,ann=None):
			
		"""Draws line into axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of domain.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],proj=['3d'])
			ax=axes[0]
		ax.plot([self.v1.x[0],self.v2.x[0]],[self.v1.x[1],self.v2.x[1]],zs=[self.v1.x[2],self.v2.x[2]],color=color,linestyle='-')
		if ann:
			m=self.getMiddle()
			ax.text(m[0]+self.domain.annXOffset, m[1]+self.domain.annYOffset, m[2]+self.domain.annZOffset, "l"+str(self.Id), None)
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
		
	def getLastVertex(self,orientation):
		
		"""Returns last vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.v2
		elif orientation==-1:
			return self.v1
		else:
			printError("Cannot return last vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def getFirstVertex(self,orientation):
		
		"""Returns first vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.v1
		elif orientation==-1:
			return self.v2
		else:
			printError("Cannot return first vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def writeToFile(self,f):
		
		"""Writes line to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Line("+str(self.Id)+")= {" + str(self.v1.Id) + "," + str(self.v2.Id) + "};\n" )
		
		return f
	
	
class arc(edge):
	
	"""Arc class storing information from gmsh .geo cicle.
	
	Will compute ``angleOffset``, ``angle`` and ``pOffset`` on creation.
	
	.. image:: ../imgs/pyfrp_gmsh_geometry/arc.png
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain arc belongs to.
		vstart (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
		vcenter (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Center vertex.
		vend (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
		Id (int): ID of arc.
		
	"""		
	
	def __init__(self,domain,vstart,vcenter,vend,Id):
		
		edge.__init__(self,domain,Id,1)
		
		self.vcenter=vcenter
		self.vstart=vstart
		self.vend=vend
		
		self.radius=self.computeRadius()
		
		self.pOffset=self.computePOffset()
		
		self.angleOffset=self.computeAngleOffset()
		self.angle=self.computeAngle()
		
	

	def computeAngleOffset(self):
		
		"""Computes and returns offset angle of arc.
		"""
		
		self.angleOffset=getAngle(self.pOffset,self.vstart.x-self.vcenter.x)
		
		return self.angleOffset
	
	def computeAngle(self):
		
		"""Computes and returns angle of arc.
		"""
		
		self.angle=getAngle(self.vstart.x-self.vcenter.x,self.vend.x-self.vcenter.x)
		return self.angle
	
	def computePOffset(self):
		
		"""Computes and returns offset point of arc.
		"""
		
		v1n,v2nb = self.getNormVec()
		
		self.pOffset=self.radius*v2nb
		self.pOffset=self.pOffset/np.linalg.norm(self.pOffset)
		
		return self.pOffset
	
	def getNormVec(self):
		
		"""Computes and returns vectors normal to arc.
		
		Returns:
			tuple: Tuple containing:
				
				* v1n (numpy.ndarray): Normal vector to ``vstart-vcenter``.
				* v2n (numpy.ndarray): Normal vector to ``vend-vcenter``.
					
		"""
		
		v1=self.vstart.x-self.vcenter.x
		v2=self.vend.x-self.vcenter.x
		
		self.v1n = v1/np.linalg.norm(v1)
		v2n = v2/np.linalg.norm(v2)
		v2nb = v2n-np.dot(v2n,self.v1n)*self.v1n
		
		self.v2nb = v2nb/np.linalg.norm(v2nb)
		
		return self.v1n,self.v2nb
	
	def getPlotVec(self):
		
		"""Returns vectors for plotting arc.
		
		Returns:
			tuple: Tuple containing:
				
				* x (numpy.ndarray): x-array.
				* y (numpy.ndarray): y-array.
				* z (numpy.ndarray): z-array.
						
		"""
		
		self.getNormVec()
			
		if np.mod(self.angle,np.pi/2.)<0.01:
			a = np.linspace(0,self.angle,1000)
		else:
			a = np.linspace(self.angleOffset-self.angle,self.angleOffset,1000)
				
		x,y,z=self.getPointOnArc(a)
		
		return x,y,z
	
	def getPointOnArc(self,a):
		
		"""Returns point on arc at angle ``a``.
		
		Returns:
			tuple: Tuple containing:
				
				* x (float): x-coordinate.
				* y (float): y-coordinate.
				* z (float): z-coordinate.
						
		"""
			
		x = self.vcenter.x[0]+np.sin(a)*self.radius*self.v1n[0]+np.cos(a)*self.radius*self.v2nb[0]
		y = self.vcenter.x[1]+np.sin(a)*self.radius*self.v1n[1]+np.cos(a)*self.radius*self.v2nb[1]
		z = self.vcenter.x[2]+np.sin(a)*self.radius*self.v1n[2]+np.cos(a)*self.radius*self.v2nb[2]	
		
		return x,y,z

	def computeRadius(self):
		
		"""Computes and returns radius of arc.
		
		Returns:
			float: Radius of arc.
		"""
		
		self.radius=np.linalg.norm(self.vstart.x-self.vcenter.x)
		return self.radius
		
	def inArc(self,x,debug=False):
		
		"""Tells if coordinate ``x`` is on arc or not.
		
		Returns:
			bool: ``True`` if on arc, ``False`` otherwise.
		"""
		
		a=self.computeAngle(array([self.radius,0])-self.vcenter.x,x-self.vcenter.x)
				
		if np.mod(a,2*np.pi)<self.angle+self.angleOffset and self.angleOffset<=np.mod(a,2*np.pi):
			return True
		else:
			return False
	
	def getRadius(self):
		
		"""Returns radius of arc."""
		
		return self.radius
	
	def getAngle(self):
		
		"""Returns angle of arc."""
		
		return self.angle
	
	def getAngleOffset(self):
		
		"""Returns offset angle of arc."""
		
		return self.angleOffset
	
	def getVstart(self):
		
		"""Returns start vertex of arc."""
		
		return self.vstart
	
	def getVend(self):
		
		"""Returns end vertex of arc."""
		
		return self.vend
	
	def getXstart(self):
		
		"""Returns start coordinate of arc."""
		
		return self.vstart.x
	
	def getXend(self):
		
		"""Returns end coordinate of arc."""
		
		return self.vend.x
	
	def getVcenter(self):
		
		"""Returns center vertex of arc."""
		
		return self.vcenter
	
	def getXcenter(self):
		
		"""Returns center coordinate of arc."""
		
		return self.vcenter.x
	
	def draw(self,ax=None,color=None,ann=None):
		
		"""Draws arc into axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of domain.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],proj=['3d'])
			ax=axes[0]
			
		x,y,z=self.getPlotVec()
		
		ax.plot(x,y,zs=z,color=color,linestyle='-')
		
		if ann:
			x,y,z=self.getPointOnArc(self.angle/2.)
			ax.text(x+self.domain.annXOffset, y+self.domain.annYOffset, z+self.domain.annZOffset, "c"+str(self.Id), None)
			
			
			
		pyfrp_plot_module.redraw(ax)
	
	def getLastVertex(self,orientation):
		
		"""Returns last vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.getVend()
		elif orientation==-1:
			return self.getVstart()
		else:
			printError("Cannot return last vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def getFirstVertex(self,orientation):
		
		"""Returns first vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==-1:
			return self.getVend()
		elif orientation==1:
			return self.getVstart()
		else:
			printError("Cannot return first vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def writeToFile(self,f):
		
		"""Writes arc to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Circle("+str(self.Id)+")= {" + str(self.vstart.Id) + ","+ str(self.vcenter.Id)+ "," + str(self.vend.Id) + "};\n" )
		
		return f
	
class lineLoop:
	
	"""Lineloop class storing information from gmsh .geo.

	Object has two major attributes:
	
		* edges (list): List of pyfrp.moduels.pyfrp_gmsh_geometry.edge objects.
		* orientations (list): List of orientations of each element, either ``1`` or ``-1`` 
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain loop belongs to.
		edgeIDs (list): List of edge IDs.
		Id (int): ID of loop.
		
	"""		
	
	def __init__(self,domain,edgeIDs,ID):
		
		self.domain=domain
		self.Id=ID
		
		self.edges,self.orientations=self.initEdges(edgeIDs)
		
	def initEdges(self,IDs):
		
		"""Constructs ``edges`` and ``orientations`` list at object initiations 
		from list of IDs.
		
		Args:
			IDs (list): List of IDs
			
		Returns:
			tuple: Tuple containing:
			
				* edges (list): List of pyfrp.moduels.pyfrp_gmsh_geometry.edge objects.
				* orientations (list): List of orientations of each element, either ``1`` or ``-1`` 
		
		"""
		
		self.edges=[]
		self.orientations=[]
		
		for ID in IDs:
			self.addEdgeByID(ID)
		
		return self.edges,self.orientations
		
	def addEdgeByID(self,ID):
		
		"""Adds edge to lineloop.
		
		Args:
			ID (int): ID of edge to be added.
			
		Returns:
			list: Updated edgeIDs list.
			
		"""
		
		
		self.edges.append(self.domain.getEdgeById(abs(ID))[0])
		self.orientations.append(np.sign(ID))
		
		return self.edges
	
	def insertEdgeByID(self,ID,pos):
		
		"""Inserts edge to lineloop at position.
		
		Args:
			ID (int): ID of edge to be inserted.
			pos (int): Position at which ID to be inserted.
			
		Returns:
			list: Updated edgeIDs list.
			
		"""
		
		self.edges.insert(pos,self.domain.getEdgeById(abs(ID))[0])
		self.orientations.insert(pos,np.sign(ID))
		
		return self.edges
	
	def removeEdgeByID(self,ID):
		
		"""Remove edge from lineloop.
		
		Args:
			ID (int): ID of edge to be removed.
			
		Returns:
			list: Updated edgeIDs list.
			
		"""
		
		idx=self.edges.index(abs(ID))
		self.edges.remove(abs(ID))
		self.orientations.pop(idx)
		
		return self.edges
	
	def reverseEdge(self,ID):
		
		"""Reverses the orientation of an edge in the line loop.
		
		Args:
			ID (int): ID of edge to be reversed.
		
		Returns:
			list: Updated orientations list.
			
		"""
		
		self.orientations[self.edges.index(abs(ID))]=-self.orientations[self.edges.index(abs(ID))]
		
		return self.orientations
		
	def checkClosed(self,debug=False):
		
		"""Checks if lineLoop is closed.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			bool: True if closed.
		
		"""
		
		b=True
		
		for i in range(len(self.edges)):
			
			#Get ID of edge
			edge1Temp=self.edges[i]
			orient1=self.orientations[i]
			
			#Get ID of next edge
			if i==len(self.edges)-1:
				edge2Temp=self.edges[i+1]
				orient2=self.orientations[i+1]
			else:
				edge2Temp=self.edges[0]
				orient2=self.orientations[0]
				
			#Get ID of first/last vertex
			firstVertexId=edge1Temp.getFirstVertex(orient1).Id
			lastVertexId=edge2Temp.getLastVertex(orient2).Id
			
			#Check if vertices are matching
			if firstVertexId!=lastVertexId:
				b=False
				
				if debug:
					printWarning("lineLoop with ID " + str(self.Id) + " does not close." )
					print "Edge with ID " +str(edge1Temp.Id) + " is not matching edge with ID " + str(edge2Temp.Id)
						
		return b
	
	def writeToFile(self,f):
		
		"""Writes line loop to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Line Loop("+str(self.Id)+")= {" )
			
		for i,s in enumerate(self.edges):
			f.write(str(self.orientations[i]*s.Id))
			if i!=len(self.edges)-1:
				f.write(",")
			else:
				f.write("};\n")
	
		return f
	
class ruledSurface:
	
	"""ruledSurface class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		loopID (int): ID of surrounding loop.
		Id (int): ID of surface.
		
	"""		
	
	def __init__(self,domain,loopID,ID):
		
		self.domain=domain
		self.Id=ID
		
		self.lineLoop=self.domain.getLineLoopById(loopID)[0]
	
	def writeToFile(self,f):
		
		"""Writes ruled surface to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Ruled Surface("+str(self.Id)+")= {"+str(self.lineLoop.Id)+ "};\n" )
	
		return f
	
class surfaceLoop:
	
	"""surfaceLoop class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain loop belongs to.
		surfaceIDs (list): List of surfaces.
		Id (int): ID of loop.
		
	"""		
	
	def __init__(self,domain,surfaceIDs,ID):
		
		self.domain=domain
		self.Id=ID
		
		self.surfaces=self.initSurfaces(surfaceIDs)	
		
	def initSurfaces(self,IDs):
		
		"""Constructs ``surfaces`` list at object initiations 
		from list of IDs.
		
		Args:
			IDs (list): List of IDs.
			
		Returns:
			list: List of pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface objects.
		
		"""
		
		self.surfaces=[]
		
		for ID in IDs:
			self.addSurfaceByID(ID)
		
		return self.surfaces
	
	def addSurfaceByID(self,ID):
		
		"""Adds surface to surfaceloop.
		
		Args:
			ID (int): ID of surface to be added.
			
		Returns:
			list: Updated surfaceIDs list.
			
		"""
		
		self.surfaces.append(self.domain.getRuledSurfaceById(ID)[0])
		
		return self.surfaces
	
	def insertSurfaceByID(self,ID,pos):
		
		"""Inserts surface to surfaceloop at position.
		
		Args:
			ID (int): ID of surface to be inserted.
			pos (int): Position at which ID to be inserted.
			
		Returns:
			list: Updated surfaceIDs list.
			
		"""
		
		self.surfaces.insert(pos,self.domain.getRuledSurfaceById(ID)[0])
		
		return self.surfaces
	
	def removeSurfaceByID(self,ID):
		
		"""Remove surface from surfaceloop.
		
		Args:
			ID (int): ID of surface to be removed.
			
		Returns:
			list: Updated surfaceIDs list.
			
		"""
		
		self.surfaces.remove(self.domain.getRuledSurfaceById(ID)[0])
		
		return self.surfaces
	
	def writeToFile(self,f):
		
		"""Writes surface loop to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Surface Loop("+str(self.Id)+")= {" )
		
		for i,s in enumerate(self.surfaces):
			f.write(str(s.Id))
			if i!=len(self.surfaces)-1:
				f.write(",")
			else:
				f.write("};\n")
	
		return f
	
class volume:

	"""Volume class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		surfaceLoopID (int): ID of surrounding surface loop.
		Id (int): ID of surface loop.
		
	"""		
	
	def __init__(self,domain,surfaceLoopID,ID):
		
		self.domain=domain
		self.Id=ID
		
		self.surfaceLoop=self.domain.getSurfaceLoopById(surfaceLoopID)[0]
	
	def writeToFile(self,f):
		
		"""Writes Volume to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Volume("+str(self.Id)+")= {"+str(self.surfaceLoop.Id)+ "};\n" )
	
		return f
				
			