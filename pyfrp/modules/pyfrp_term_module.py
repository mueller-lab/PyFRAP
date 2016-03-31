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

import colorama
import PyQt4.QtGui as QtGui

def printWarning(txt):

	print(colorama.Fore.YELLOW + "WARNING:") + colorama.Fore.RESET + txt

def printError(txt):
	
	print(colorama.Fore.RED + "ERROR:") + colorama.Fore.RESET + txt

def printNote(txt):

	print(colorama.Fore.GREEN + "NOTE:") + colorama.Fore.RESET + txt
	
def printDict(dic):
	for k in dic.keys():
		print k , ' = ' , dic[k]
	return True	

def printObjAttr(var,obj):
        print var, " = ", vars(obj)[str(var)]
        return var

def printAllObjAttr(obj):
	for item in vars(obj):
		print item, " = ", vars(obj)[str(item)]
	return True	