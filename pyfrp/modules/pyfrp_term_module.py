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

"""Terminal module for PyFRAP toolbox. Provides extra functions for a nicer 
custom output inside a Python/bash terminal.
"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

import colorama
import PyQt4.QtGui as QtGui
import numpy as np

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================


def printWarning(txt):
	
	"""Prints Warning of the form "WARNING: txt", while warning is rendered yellow.	
	"""

	print(colorama.Fore.YELLOW + "WARNING:") + colorama.Fore.RESET + txt

def printError(txt):
	
	"""Prints Error of the form "ERROR: txt", while error is rendered red.	
	"""
	
	print(colorama.Fore.RED + "ERROR:") + colorama.Fore.RESET + txt

def printNote(txt):
	
	"""Prints note of the form "NOTE: txt", while note is rendered green.	
	"""

	print(colorama.Fore.GREEN + "NOTE:") + colorama.Fore.RESET + txt
	
def printDict(dic,maxL=5):
	
	"""Prints all dictionary entries in the form key = value.
	
	If attributes are of type ``list`` or ``numpy.ndarray``, will check if the size
	exceeds threshhold. If so, will only print type and dimension of attribute.
	
	Args:
		dic (dict): Dictionary to be printed.
		
	Returns:
		bool: True
	
	"""
	
	for k in dic.keys():
		printAttr(k,dic[k],maxL=maxL)
		
	return True	

def printObjAttr(var,obj):
	
	"""Prints single object attribute in the form attributeName = attributeValue.
	Args:
		var (str): Name of attribute.
		obj (object): Object to be printed.
	
	Returns:
		str: Name of attribute.
	
	"""
	
	
        print var, " = ", vars(obj)[str(var)]
        return var

def printAllObjAttr(obj,maxL=5):
	
	"""Prints all object attributes in the form attributeName = attributeValue.
	
	If attributes are of type ``list`` or ``numpy.ndarray``, will check if the size
	exceeds threshhold. If so, will only print type and dimension of attribute.
	
	Args:
		obj (object): Object to be printed.
	
	Keyword Args:
		maxL (int): Maximum length threshhold.
	
	"""
	
	for item in vars(obj):
		printAttr(item,vars(obj)[str(item)],maxL=maxL)
		
	return True

def printAttr(name,attr,maxL=5):

	"""Prints single attribute in the form attributeName = attributeValue.
	
	If attributes are of type ``list`` or ``numpy.ndarray``, will check if the size
	exceeds threshhold. If so, will only print type and dimension of attribute.
	
	Args:
		name (str): Name of attribute.
		attr (any): Attribute value.
		
	Keyword Args:
		maxL (int): Maximum length threshhold.
	
	"""

	if isinstance(attr,(list)):
		if len(attr)>maxL:
			print name, " = ", getListDetailsString(attr)
			return True
	elif isinstance(attr,(np.ndarray)):
		if min(attr.shape)>maxL:
			print name, " = ", getArrayDetailsString(attr)
			return True
		
	print name, " = ", attr
		
	return True	


def getListDetailsString(l):
	
	"""Returns string saying "List of length x", where x is the length of the list. 
	
	Args:
		l (list): Some list.
	
	Returns:
		str: Printout of type and length.
	"""
	
	return "List of length " + str(len(l))

def getArrayDetailsString(l):
		
	"""Returns string saying "Array of shape x", where x is the shape of the array. 
	
	Args:
		l (numpy.ndarray): Some array.
	
	Returns:
		str: Printout of type and shape.
	"""
	
	return "Array of shape " + str(l.shape)

