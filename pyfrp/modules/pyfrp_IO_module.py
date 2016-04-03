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

#Input/Output module for PyFRAP toolbox, including following functions:

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

import pickle
import platform
import gc
import sys
import os

from pyfrp.modules import pyfrp_misc_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def saveToPickle(obj,fn=None):
	
	"""Saves obj into pickled format.
	
	.. note:: If ``fn==Non``, will try to save to ``obj.name``, otherwise unnamed.pk
	
	Keyword Args:
		fn (str): Output file name.	
	
	Returns: 
		str: Output filename.
	
	"""
	
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
	
	"""Loads obj from pickled format.
	
	Args:
		fn (str): Filename.	
	
	Returns: 
		str: Output filename.
	
	"""
	
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
	
	"""Loads molecule object from pickle file
	and brings it up-to-date.
	
	Args:
		fn (str): Filename.	
	
	Keyword Args: 
		update (bool): Update to current version.
	
	Returns: 
		pyfrp.subclasses.pyfrp_molecule: Molecule file.
	
	"""
	
	mol=loadFromPickle(fn)
	if update:
		mol.update_version()
	return mol

def loadEmbryo(fn,update=True):
	
	"""Loads embryo object from pickle file
	and brings it up-to-date.
	
	Args:
		fn (str): Filename.	
	
	Keyword Args: 
		update (bool): Update to current version.
	
	Returns: 
		pyfrp.subclasses.pyfrp_embryo: Embryo file.
	
	"""
	
	emb=loadFromPickle(fn)
	if update:
		emb.update_version()
	return emb

def cleanUp():
	"""Calls garbage collector to clean up.
	"""
	
	gc.collect()
	return None


