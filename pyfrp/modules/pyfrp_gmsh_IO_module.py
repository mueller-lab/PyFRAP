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

"""PyFRAP module for reading/writing gmsh .geo files. Module mainly has the following features:

	* Read .geo files.
	* Translate geometric entities and variables defined in .geo files.
	* Construct :py:class:`pyfrp.pyfrp_gmsh_geometry.domain` object describing complete geometry.
	* Update parameters in .geo files.
	* Add/Remove some geometric entities.
	* Add/update box fields to allow refinement of certain ROIs in mesh.

This module together with pyfrp.pyfrp_gmsh_geometry and pyfrp.pyfrp_gmsh_module works partially as a python gmsh wrapper, however is incomplete.
If you want to know more about gmsh, go to http://gmsh.info/doc/texinfo/gmsh.html .
	
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy (use indirect import here, so convertMathExpr relates to numpy functions automatically when translating)
from numpy import *

#PyFRAP Modules
import pyfrp_gmsh_geometry
import pyfrp_misc_module

#Misc
import shutil
from tempfile import mkstemp
import os

                   
#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================


def splitLine(line,delim="=",closer=";"):
	
	"""Splits line at ``delim``, trimming ``closer``.
	
	Example:
	
	>>> splitLine("Point(3)={1,3,1};")
	>>> ("Point(3)","{1,3,1}")
	
	Args:
		line (str): Line to be splitted.
	
	Keyword Args: 
		delim (str): Delimiter at which to be splitted.
		closer (str): Closing character to be trimmed.
		
	Returns:
		tuple: Tuple containing:
			
			* var (str): Name of variable.
			* val (str): Value of variable
		
	"""
	
	line=line.strip()
	var,val=line.split(delim)
	
	val=val.strip('\n')
	val=val.strip(closer)
	
	var=var.strip()
	val=val.strip()
	
	return var,val

def getId(var,delimOpen="(",delimClose=")"):
	
	"""Returns ID of object that is given between the delimiters ``delimOpen``
	and ``delimClose``.
	
	Example:
	
	>>> getId("Point(3)")
	>>> ("Point",3) 
	
	Args:
		var (str): String describing geoFile variable name.
	
	Keyword Args: 
		delimOpen (str): Openening delimiter of ID.
		delimClose (str): Closing delimiter of ID.
		
	Returns:
		tuple: Tuple containing:
			
			* typ (str): Type of geometric variable.
			* Id (str): Id of geometric variable.
		
	"""
	
	typ,Id=var.split(delimOpen)
	Id=int(Id.split(delimClose)[0])
	return typ,Id

def getVals(val,parmDic):
	
	"""Translates value of parameter into list of floats.
	
	Uses parameter dictionary to translate predefined variables into
	floats.
	
	Example:
	
	>>> getVals("{10,3,5}")
	>>> [10.,3.,5.]
		
	Args:
		val (str): Value string of geometric variable.
		parmDic (dict): Parameter dictionary.
	
	Returns:
		rList (list): List of translated values.
	
	"""
	
	valList,i=pyfrp_misc_module.str2list(val,dtype="str",openDelim="{",closeDelim="}",sep=",")
	rList=[]
	
	for v in valList:
		v=v.strip()
		rList.append(applyParmDic(v,parmDic))
	return rList

def convertMathExpr(val):
	
	"""Converts math expressions from .geo syntax into python
	syntax.
	
	.. note::Not all translations have been implemented yet. You can 
	   simply add here expressions by adding a translation to the 
	   translations list (``translations.append([CExpression,PythonExpression])``).
	
	"""
	
	translations=[]
	translations.append(["^","**"])
	translations.append(["Sqrt","sqrt"])
	
	for translation in translations:
		val=val.replace(translation[0],translation[1])
	
	return val

def applyParmDic(val,parmDic):
	
	"""Applies parameter dictionary to variable value.
	
	Example: 
	
	>>> parmDic={'radius',3}
	>>> applyParmDic('radius',parmDic)
	>>> 3
	
	And also applies mathemtical expressions:
	
	>>> parmDic={'radius',3}
	>>> applyParmDic('radius^2-radius',parmDic)
	>>> 6
	
	Args:
		val (str): Value string of geometric variable.
		parmDic (dict): Parameter dictionary.
		
	Returns:
		val (float): Evaluated value.
	
	"""
	
	keys,lengths=sortKeysByLength(parmDic)
	keys.reverse()
	
	val=convertMathExpr(val)
	
	for key in keys:
		val=val.replace(key,str(parmDic[key]))	
	
	return eval(val)
	
def sortKeysByLength(dic):
	
	"""Sorts dictionary by length of keys. 
	"""
	
	lengths=[]
	for key in dic.keys():
		lengths.append(len(key))
	
	keys,lengths=pyfrp_misc_module.sortListsWithKey(dic.keys(),lengths)
	return keys,lengths
		
def readParameter(line,parmDic):
	
	"""Reads in parameter from line and translates values using ``parmDic``.
	
	Args:
		line (str): Line to be splitted.
		parmDic (dict): Parameter dictionary.
		
	Returns:
		tuple: Tuple containing:
		
			* var (str): Name of variable.
			* val (float): Value of variable
	
	"""
	
	var,val = splitLine(line)
	val=applyParmDic(val,parmDic)
	return var,val

def readLine(line,parmDic,domain):
	
	"""Reads in line from .geo file. 
	
	Tries to extract type of geometric object and its parameters 
	and uses this to append a geomtric entity to ``domain``. 
	
	If ``line`` describes a parameter, stores parameter name and its value
	in ``parmDic``.
	
	Args:
		line (str): Line to be splitted.
		parmDic (dict): Parameter dictionary.
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object, storing all geometric entities.
		
	Returns:
		tuple: Tuple containing:
		
			* parmDic (dict): Updated parameter dictionary.
			* typ (str): Object type type.
			* Id (int): ID of object.
			* vals (list): Values of object.
			* domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Updated domain object.
			
	"""
	
	if line.startswith('//'):
		#This line is a comment, return parmDic
		return parmDic, "comment", -1, [],domain
	
	if "{" in line or "}" in line:
		#This line is some sort of object
		
		var,val = splitLine(line)
		typ,Id = getId(var)
		vals=getVals(val,parmDic)
		
		if typ=="Point":
			domain.addVertex([vals[0],vals[1],vals[2]],Id=Id,volSize=vals[3])
		elif typ=="Line":
			
			v1,idx=domain.getVertexById(vals[0])
			v2,idx=domain.getVertexById(vals[1])
			domain.addLine(v1,v2,Id=Id)
		elif typ=="Circle":
			vstart,idx=domain.getVertexById(vals[0])		
			vcenter,idx=domain.getVertexById(vals[1])
			vend,idx=domain.getVertexById(vals[2])
			domain.addArc(vstart,vcenter,vend,Id=Id)
		else:
			#This is a object like surface or volume which we don't care about
			pass
		
		return parmDic, typ, Id , vals,domain
			
	else:
		#Check if line is empty
		if len(line.strip())>0:
			
			#This is a parameter line, append new parameter to parmDic
			var,val=readParameter(line,parmDic)
			
			parmDic[var]=val
			return parmDic, "parameter", -2, [],domain
		
		else:
			return parmDic, "empty", -3, [], domain
		
def readGeoFile(fn):
		
	"""Reads in .geo file and tries to extract geometry defined in .geo file
	into a :py:class`pyfrp.modules.pyfrp_gmsh_geometry.domain`.
	
	Args:
		fn (str): Filename of .geo file.
	
	Returns:
		tuple: Tuple containing:
		
			* parmDic (dict): Updated parameter dictionary.
			* domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object.
			
	"""
	
	#new parameter dictionary
	parmDic={}
	
	#New domain
	domain=pyfrp_gmsh_geometry.domain()
	
	#Read file
	with open(fn,'r') as f:
		for line in f:
			parmDic,typ,Id,vals,domain=readLine(line,parmDic,domain)
			
	return domain,parmDic
			
def txtLineReplace(filePath, pattern, subst):
		
	"""Replaces line in file that starts with ``pattern`` and substitutes it 
	with ``subst``.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	Args:
		filePath (str): Filename.
		pattern (str): Pattern to be looked for.
		subst (str): String used as a replacement.
			
	"""
	
	
	#Create temp file
	fh, absPath = mkstemp()
	newFile = open(absPath,'w')
	oldFile = open(filePath)
	
	#Loop through file and replace line 
	for line in oldFile:
		
		if line.startswith(pattern):
			newFile.write(line.replace(line, subst))
		else:
			newFile.write(line)
			
	#close temp file
	newFile.close()
	os.close(fh)
	oldFile.close()
		
	#Remove original file
	os.remove(filePath)
	
	#Move new file
	shutil.move(absPath, filePath)
	return			

def updateParmGeoFile(fn,name,val):
	
	"""Updates parameter in .geo file.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	Args:
		fn (str): Filename of .geo file.
		name (str): Name of parameter.
		val (float): Value of parameter.
			
	"""
		
	
	substr=name+"="+str(val)+";"+'\n'
	txtLineReplace(fn,name,substr)
	
	return

def getAllIDsOfType(fn,elementType):
	
	"""Finds all IDs of a specific .geo element type in a .geo file.
	
	Args:
		fn (str): Filename of .geo file.
		elementType (str): Type of parameter, for example ``"Point"``.
		
	Returns:
		list: List of IDs.
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return
	
	Ids=[]
	
	f = open (fn,'rb')
	for line in f:
		
		if line.strip().startswith(elementType):
			
			var,val = splitLine(line)
			if elementType=="Field":
				typ,Id = getId(var,delimOpen="[",delimClose="]")
			else:
				typ,Id = getId(var)
			
			Ids.append(Id)
	Ids=list(unique(Ids))			
	
	f.close()
	return Ids
	
def getLargestIDOfType(fn,elementType):
	
	"""Finds largest ID of a specific .geo element type in a .geo file.
	
	Args:
		fn (str): Filename of .geo file.
		elementType (str): Type of parameter, for example ``"Point"``.
		
	Returns:
		int: Largest ID.
	"""
	
	return max(getAllIDsOfType(fn,elementType))

def getBkgdFieldID(fn):
	
	"""Finds ID of background field in .geo file.
	
	.. note:: Will return ``None`` if .geo file has no background
	   field specified.
	
	Args:
		fn (str): Filename of .geo file.
			
	Returns:
		int: ID of background field.
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return
	
	f = open (fn,'rb')
	for line in f:
		if line.strip().startswith("Background Field"):
			var,val = splitLine(line)
			f.close()
			return int(val.strip())
	
	
	f.close()	
	return 

def getLastNonEmptyLine(fn):
	
	"""Finds index of last non-empty line in .geo file.
	
	Args:
		fn (str): Filename of .geo file.
			
	Returns:
		int: Index of last non-empty line.
	"""
	
	idx=0
	with open(fn,'rb') as f:
		for i,line in enumerate(f):
			if len(line.strip()):
				idx=i
	return idx
	
def removeTailingLines(filePath,idx):
	
	"""Removes all empty lines at the end of a .geo file.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	Args:
		filePath (str): Filename of .geo file.
		idx (int): Index of last non-empty line
			
	"""
	
	#Create temp file
	fh, absPath = mkstemp()
	newFile = open(absPath,'w')
	oldFile = open(filePath)
	
	#Loop through file and only write until line idx
	for i,line in enumerate(oldFile):
		
		if i<=idx:
			newFile.write(line)
			
	#close temp file
	newFile.close()
	os.close(fh)
	oldFile.close()
		
	#Remove original file
	os.remove(filePath)
	
	#Move new file
	shutil.move(absPath, filePath)
	return		
	
def copyIntoTempFile(fn,close=True):
	
	"""Copies file into tempfile.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	.. note:: If ``close==True``, will return ``fh=None`` and ``tempFile=None``.
	
	Args:
		fn (str): Filename of file.
	
	Keyword Args:
		close (bool): Close files after copying.
	
	Returns:
		tuple: Tuple containing:
		
			* tempFile (file): File handle to temp file.
			* fh (file): File handle to original file.
			* tempPath (tempPath): Path to temp file.
	
	"""
	
	oldFile = open(fn)
	fh, tempPath = mkstemp()
	tempFile = open(tempPath,'w')
	
	for line in oldFile:
		tempFile.write(line)

	oldFile.close()
	
	if close:
		tempFile.close()
		os.close(fh)
		tempFile=None
		fh=None
		
	return tempFile, fh,tempPath

def getLinesByID(fn,elementId,elementType=""):
	
	"""Finds all lines in .geo file that contain geometric entitity with ID ``elementId``.
	
	.. note:: IDs in geometric files can be given per entitity type. That is, one can have
	   for example a point with ID=1 (``Point(1)``) aswell as a line with ID=1 (``Line(1)``).
	   Thus one may want to use ``elementType`` to restrict the search for a specific element type.
	
	Args:
		fn (str): Filename of .geo file.
		elementId (int): ID to look for.
	
	Keyword Args:
		elementType (str): Type of element to restrict search on.
	
	Returns:
		list: Line numbers at which element appears.
	
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return
	
	f = open (fn,'rb')
	
	lineNumbers=[]
	
	for i,line in enumerate(f):
		if line.strip().startswith(elementType):
			
			var,val = splitLine(line)
			if elementType=="Field":
				typ,Id = getId(var,delimOpen="[",delimClose="]")
			else:
				typ,Id = getId(var)
			
			if elementId==Id:
				lineNumbers.append(i)
	
	f.close()
	return lineNumbers
				
def removeElementFromFile(fn,elementType,elementId,delimOpen="(",delimClose=")"):
	
	"""Removes element with type ``elementType`` and ID ``elementID`` from .geo file.
	
	Args:
		fn (str): Filename of .geo file.
		elementId (int): ID of element to remove.
		elementType (str): Type of element to remove.
	
	Keyword Args:
		delimOpen (str): Openening delimiter of ID.
		delimClose (str): Closing delimiter of ID.
	

	"""
	
	txtLineReplace(fn,elementType+delimOpen+str(elementId)+delimClose,"")
	return

def addBoxField(fn,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ,comment="",fnOut=""):		
	
	"""Adds box field to .geo file by doing the following:
		
		* Copies file into temp file for backup using :py:func:`copyIntoTempFile`.
		* Finds all IDs of previous ``Field`` entities using :py:func:`getAllIDsOfType`.
		* Finds current background field using :py:func:`getBkgdFieldID` .
		* If previous fields exist, removes them from file using :py:func:`removeElementFromFile` .
		* Finds last non-empty line using  :py:func:`getLastNonEmptyLine` .
		* Removes empty lines at end of file using :py:func:`removeTailingLines` .
		* Writes comment using  :py:func:`writeComment` .
		* Writes box field using  :py:func:`writeBoxField` .
		* Writes background field using  :py:func:`writeBackgroundField` .
		
	.. note:: Comment is useful to describe in .geo file what the the box field actually does.
	
	.. note:: Generally, background field will use ``volSizeIn`` as background mesh volume size.
	
	.. note:: Unit for parameter is pixels.
	
	.. note:: If ``fnOut`` is not specified, will overwrite input file.
	
	See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		fn (str): Filename of .geo file.
		volSizeIn (float): Mesh element volume inside box.
		volSizeOut (float): Mesh element volume outside box.
		rangeX (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
		rangeY (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
		rangeZ (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
		
	Keyword Args:
		comment (str): Comment to be added before box field.
		fnOut (str): Filepath for output.
	

	"""
	
	#Copy everything into tempfile
	tempFile,fh,tempPath = copyIntoTempFile(fn,close=True)
	
	#Find if there are already Fields defined
	fieldIDs=getAllIDsOfType(tempPath,"Field")
	
	#Find background fields
	bkgdID=getBkgdFieldID(tempPath)
	
	#Remove bkgdID from field IDs
	if bkgdID!=None:
		fieldIDs.remove(bkgdID)
	
	#If there is already a background field, remove all lines containing it
	if bkgdID!=None:
		removeElementFromFile(tempPath,"Field",bkgdID,delimOpen="[",delimClose="]")
		removeElementFromFile(tempPath,"Mesh.",bkgdID,delimOpen="",delimClose="")
		removeElementFromFile(tempPath,"Background Field =",bkgdID,delimOpen="",delimClose="")
	
	#Get index of last non-empty line
	idxNonEmpty=getLastNonEmptyLine(tempPath)
	
	#remove tailing lines
	removeTailingLines(tempPath,idxNonEmpty)
	
	#Open file again for appending new lines
	with open(tempPath,'a') as f:
		
		#Write Empty line
		f.write('\n')
		
		#Write Comment
		f=writeComment(f,comment)
		
		#Get ID of new field
		if len(fieldIDs)>0:
			newFieldID=max(fieldIDs)+1
		else:
			newFieldID=1
		
		#Write new box field
		f=writeBoxField(f,newFieldID,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ)
		
		#Append to field ids
		fieldIDs.append(newFieldID)
		
		#Write Background field
		f=writeBackgroundField(f,max(fieldIDs)+1,fieldIDs)
	
	#Move new file either to fnOut or to fn itself
	if fnOut!="":
		shutil.move(tempPath, fnOut)
	else:
		shutil.move(tempPath, fn)
	
	return	
	
def writeComment(f,comment):
	
	"""Writes comment line into file.
	
	Args:
		f (file): Filehandle.
		comment (str): Comment to be written.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("//"+comment+"\n")
	return f
	
def writeBackgroundField(f,fieldID,ids):
	
	"""Writes background field into into file.
	
	.. note:: Will take finest mesh for background field.
	   See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new background field.
		ids (list): List of field IDs used for background mesh computation.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;"+"\n")
	f.write("Field["+str(fieldID)+"] = Min"+";"+"\n")
	f.write("Field["+str(fieldID)+"].FieldsList = {")
	for i,d in enumerate(ids):
		f.write(str(d))
		if i<len(ids)-1:
			f.write(",")
	f.write("}"+";"+"\n")
	f.write('\n')
	f.write("Background Field ="+str(fieldID) +";"+"\n")
	f.write('\n')
	return f
	
	
def writeBoxField(f,fieldID,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ):	
	
	"""Writes box field into into file.
	
	See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new box field.
		volSizeIn (float): Mesh element volume inside box.
		volSizeOut (float): Mesh element volume outside box.
		rangeX (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
		rangeY (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
		rangeZ (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
		
	Returns:
		file: Filehandle.
	
	"""
	
	
	f.write("Field["+str(fieldID)+"] = Box"+";"+"\n")
	f.write("Field["+str(fieldID)+"].VIn = "+str(volSizeIn)+";"+"\n")
	f.write("Field["+str(fieldID)+"].VOut = "+str(volSizeOut)+";"+"\n")
	f.write("Field["+str(fieldID)+"].XMin = "+str(rangeX[0])+";"+"\n")
	f.write("Field["+str(fieldID)+"].XMax = "+str(rangeX[1])+";"+"\n")
	f.write("Field["+str(fieldID)+"].YMin = "+str(rangeY[0])+";"+"\n")
	f.write("Field["+str(fieldID)+"].YMax = "+str(rangeY[1])+";"+"\n")
	f.write("Field["+str(fieldID)+"].ZMin = "+str(rangeZ[0])+";"+"\n")
	f.write("Field["+str(fieldID)+"].ZMax = "+str(rangeZ[1])+";"+"\n")
	f.write('\n')
	return f
	
	


		

	