#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Input/Output module for PyFRAP toolbox, including following functions:

#


#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

from numpy import *
import csv
import time
import os
import pickle
import platform

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def saveToPickle(obj,fn=None):

        if fn==None:
                if hasattr(obj,"name"):
                        fn=obj.name+".pk"
                else:
                        fn="unnamed"+".pk"
                
        with open(fn, 'wb') as output:
                pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        
        return fn

def loadFromPickle(fn):
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



