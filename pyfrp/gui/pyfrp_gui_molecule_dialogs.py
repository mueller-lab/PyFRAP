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

#PyQT Dialogs for editing molecule objects
#(1) moleculeDialog

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#QT
from PyQt4 import QtGui, QtCore

#===================================================================================================================================
#Dialog for select/edit molecule
#===================================================================================================================================

class moleculeDialog(QtGui.QDialog):
	
	def __init__(self,molecule,parent):
		
		super(moleculeDialog,self).__init__(parent)
		
		#Pass molecule
		self.molecule=molecule
		
		#Buttons
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Labels
		self.lblName = QtGui.QLabel("Name:", self)
		
		#LineEdit
		self.qleName = QtGui.QLineEdit(self.molecule.name)
		self.qleName.editingFinished.connect(self.setName)
		
		#Layout
		grid = QtGui.QGridLayout()
		grid.addWidget(self.lblName,1,1)
		grid.addWidget(self.qleName,2,1)
		grid.addWidget(self.btnDone,2,2)
		
		self.setLayout(grid)    
			
		self.setWindowTitle('Edit Molecule')    
		self.show()
		
	def setName(self):
		text=self.qleName.text()
		self.molecule.setName(str(text))
		
		return self.molecule.getName()
	
	def donePressed(self):
		self.done(1)
		return self.molecule