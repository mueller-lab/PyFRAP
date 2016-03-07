#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Input/Output module for PyFRAP toolbox, including following functions:




#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

import pickle
import platform
import gc
import sys
import inspect
import os

import pyfrp_misc_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def saveToPickle(obj,fn=None):
	cleanUp()
        if fn==None:
                if hasattr(obj,"name"):
                        fn=obj.name+".pk"
                else:
                        fn="unnamed"+".pk"
                
        with open(fn, 'wb') as output:
                pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        
        return fn

def loadFromPickle(fn):
	
	cleanUp()
	
	#Need to do append subclasses folder here. Sometimes pickle has problem finding the classes
	
	sys.path.append(pyfrp_misc_module.getSubclassesDir()+'/')

        if platform.system() in ["Darwin","Linux"]:
                filehandler=open(fn, 'r')
        elif platform.system() in ["Windows"]:
                filehandler=open(fn, 'rb')
                
        loadedFile=pickle.load(filehandler)
        
        return loadedFile

def loadMolecule(fn,update=True):
	mol=loadFromPickle(fn)
	if update:
		mol.update_version()
	return mol

def loadEmbryo(fn,update=True):
	emb=loadFromPickle(fn)
	if update:
		emb.update_version()
	return emb

def cleanUp():
	gc.collect()
	return None


