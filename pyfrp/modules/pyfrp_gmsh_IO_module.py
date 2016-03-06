#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Gmsh module for PyFRAP toolbox, including following functions:

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


def splitLine(line,delim="=",closer=";"):
	line=line.strip()
	var,val=line.split(delim)
	
	val=val.strip('\n')
	val=val.strip(closer)
	
	var=var.strip()
	val=val.strip()
	
	return var,val

def getId(var,delimOpen="(",delimClose=")"):
	typ,Id=var.split(delimOpen)
	Id=int(Id.split(delimClose)[0])
	return typ,Id

def getVals(val,parmDic):
	valList,i=pyfrp_misc_module.str2list(val,dtype="str",openDelim="{",closeDelim="}",sep=",")
	rList=[]
	
	for v in valList:
		v=v.strip()
		rList.append(applyParmDic(v,parmDic))
	return rList

def convertMathExpr(val):
	val=val.replace("^","**")
	val=val.replace("Sqrt","sqrt")
	
	return val

def applyParmDic(val,parmDic):
	
	keys,lengths=sortKeysByLength(parmDic)
	keys.reverse()
	
	val=convertMathExpr(val)
	
	for key in keys:
		
		val=val.replace(key,str(parmDic[key]))	
	
	return eval(val)
	
def sortKeysByLength(dic):
	lengths=[]
	for key in dic.keys():
		lengths.append(len(key))
	
	keys,lengths=pyfrp_misc_module.sortListsWithKey(dic.keys(),lengths)
	return keys,lengths
		
def readParameter(line,parmDic):
	var,val = splitLine(line)
	val=applyParmDic(val,parmDic)
	return var,val

def readLine(line,parmDic,domain):
	
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
	
	substr=name+"="+str(val)+";"+'\n'
	txtLineReplace(fn,name,substr)
	
	return

def getAllIDsOfType(fn,elementType):
	
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
	return max(getAllIDsOfType(fn,elementType))

def getBkgdFieldID(fn):
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
	idx=0
	with open(fn,'rb') as f:
		for i,line in enumerate(f):
			if len(line.strip()):
				idx=i
	return idx
	
def removeTailingLines(filePath,idx):
	
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
	txtLineReplace(fn,elementType+delimOpen+str(elementId)+delimClose,"")
	return

def addBoxField(fn,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ,comment="",fnOut=""):		
	
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
	f.write("//"+comment+"\n")
	return f
	
def writeBackgroundField(f,fieldID,ids):
	
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
	
	


		

	