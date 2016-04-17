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

#PyFRAP modules
import pyfrp_plot_module
from pyfrp_term_module import *
import pyfrp_misc_module

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

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================
	
class domain:
	
	"""Domain class storing embryo geometry entities.
	
	Args:
		edges (list): List of edges.
		vertices (list): List of edges.
		arcs (list): List of edges.
		lines (list): List of edges.
		annXOffset (float): Offset of annotations in x-direction.
		annYOffset (float): Offset of annotations in y-direction.
		annZOffset (float): Offset of annotations in z-direction.
			
	"""
	
	def __init__(self):
		
		self.edges=[]
		self.vertices=[]
		self.arcs=[]
		self.lines=[]
		
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
		
	