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

#Classes for creating/extracting gmsh geometries for PyFRAP toolbox

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
		a=np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
		
		if a<0:
			return self.getAngle(vec2,vec1)
		return a

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================
	
class domain:
	
	def __init__(self):
		
		self.edges=[]
		self.vertices=[]
		self.arcs=[]
		self.lines=[]
		
		self.annXOffset=3.
		self.annYOffset=3.
		self.annZOffset=3.
		
		
	def addVertex(self,x,Id=None,volSize=None):
		
		newId=self.getNewId(self.vertices,Id)
		
		v=vertex(self,x,newId,volSize=volSize)
		self.vertices.append(v)	
		
		return v
	
	def addLine(self,v1,v2,Id=None):
		
		newId=self.getNewId(self.lines,Id)
		
		e=line(self,v1,v2,newId)
		self.lines.append(e)
		self.edges.append(e)
		
		return e
	
	def addArc(self,vstart,vcenter,vend,Id=None):
		
		newId=self.getNewId(self.arcs,Id)
			
		a=arc(self,vstart,vcenter,vend,newId)
		self.arcs.append(a)
		self.edges.append(a)
		
		return a
	
	def checkIdExists(self,Id,objList):
		IdList=pyfrp_misc_module.objAttrToList(objList,'Id')
		if Id in IdList:
			printWarning("Object with Id " + str(Id) + " already exists.")
			return True
		return False
	
	def getNewId(self,objList,Id=None):
		
		if Id==None:
			newId=self.incrementID(objList)
		else:
			if self.checkIdExists(Id,objList):
				newId=self.incrementID(objList)
			else:
				newId=Id
		
		return newId
		
	def incrementID(self,objList):
		if len(objList)==0:
			newId=0
		else:
			IdList=pyfrp_misc_module.objAttrToList(objList,'Id')
			newId=max(IdList)+1		
		return newId
		
	def getEdgeById(self,ID):
		for i,e in enumerate(self.edges):
			if e.Id==ID:
				return e,i
		return False,False
	
	def getVertexById(self,ID):
		for i,v in enumerate(self.vertices):
			if v.Id==ID:
				return v,i
		return False,False
		
	def draw(self,ax=None,color=None,ann=None):
		
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
		
	def getExtend(self):
		x=[]
		y=[]
		z=[]
		for v in self.vertices:
			x.append(v.x[0])
			y.append(v.x[1])
			z.append(v.x[2])
		return min(x), max(x), min(y),max(y), min(z),max(z)
	
	def verticesCoordsToList(self):
		l=[]
		for v in self.vertices:
			l.append(v.x)
		return l
	
class vertex:
	
	def __init__(self,domain,x,Id,volSize=None):
		self.domain=domain
		
		self.x=np.array(x)
		self.Id=Id
		self.volSize=volSize
			
	def draw(self,ax=None,color=None,ann=None):
		
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
		self.x=x
		return self.x

class edge:
	
	def __init__(self,domain,Id,typ):
		self.domain=domain
		self.Id=Id
		self.typ=typ
		
	def getDomain(self):
		return self.domain
	
	def getId(self):
		return self.Id
	
	def getTyp(self):
		return self.typ
	
	def decodeTyp(self):
		if typ==1:
			return "arc"
		elif typ==0:
			return "line"
	
class line(edge):
	
	def __init__(self,domain,v1,v2,Id):
		
		edge.__init__(self,domain,Id,0)

		self.v1=v1
		self.v2=v2
	
	def getMiddle(self):
		return (self.v1.x+self.v2.x)/2.
	
	def draw(self,ax=None,color=None,ann=None):
		
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
	
		self.angleOffset=getAngle(self.pOffset,self.vstart.x-self.vcenter.x)
		
		return self.angleOffset
	
	def computeAngle(self):
		self.angle=getAngle(self.vstart.x-self.vcenter.x,self.vend.x-self.vcenter.x)
		return self.angle
	
	def computePOffset(self):
		
		v1n,v2nb = self.getNormVec()
		
		self.pOffset=self.radius*v2nb
		self.pOffset=self.pOffset/np.linalg.norm(self.pOffset)
		
		return self.pOffset
	
	def getNormVec(self):
		
		v1=self.vstart.x-self.vcenter.x
		v2=self.vend.x-self.vcenter.x
		
		self.v1n = v1/np.linalg.norm(v1)
		v2n = v2/np.linalg.norm(v2)
		v2nb = v2n-np.dot(v2n,self.v1n)*self.v1n
		
		self.v2nb = v2nb/np.linalg.norm(v2nb)
		
		return self.v1n,self.v2nb
	
	def getPlotVec(self):
	
		self.getNormVec()
		
	
			
		if np.mod(self.angle,np.pi/2.)<0.01:
			a = np.linspace(0,self.angle,1000)
		else:
			a = np.linspace(self.angleOffset-self.angle,self.angleOffset,1000)
			
				
		
			
		x,y,z=self.getPointOnArc(a)
		
		return x,y,z
	
	def getPointOnArc(self,a):
			
		x = self.vcenter.x[0]+np.sin(a)*self.radius*self.v1n[0]+np.cos(a)*self.radius*self.v2nb[0]
		y = self.vcenter.x[1]+np.sin(a)*self.radius*self.v1n[1]+np.cos(a)*self.radius*self.v2nb[1]
		z = self.vcenter.x[2]+np.sin(a)*self.radius*self.v1n[2]+np.cos(a)*self.radius*self.v2nb[2]	
		
		return x,y,z

	def computeRadius(self):
		self.radius=np.linalg.norm(self.vstart.x-self.vcenter.x)
		return self.radius
		
	def inArc(self,x,debug=False):
		
		a=self.computeAngle(array([self.radius,0])-self.vcenter.x,x-self.vcenter.x)
				
		if np.mod(a,2*np.pi)<self.angle+self.angleOffset and self.angleOffset<=np.mod(a,2*np.pi):
			return True
		else:
			return False
	
	def getRadius(self):
		return self.radius
	
	def getAngle(self):
		return self.angle
	
	def getAngleOffset(self):
		return self.angleOffset
	
	def getVstart(self):
		return self.vstart
	
	def getVend(self):
		return self.vend
	
	def getXstart(self):
		return self.vstart.x
	
	def getXend(self):
		return self.vend.x
	
	def getVcenter(self):
		return self.vcenter
	
	def getXcenter(self):
		return self.vcenter.x
	
	def draw(self,ax=None,color=None,ann=None):
		
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
		
	