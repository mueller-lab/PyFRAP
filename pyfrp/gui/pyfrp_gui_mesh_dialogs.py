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

#PyQT Dialogs for mesh class
#(1) 

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#QT
from PyQt4 import QtGui, QtCore

#PyFRAP GUI classes
import pyfrp_gui_basics

#PyFRAP modules
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_img_module
from pyfrp.modules import pyfrp_misc_module

#Numpy/Scipy
import numpy as np

#Misc 
import os


#===================================================================================================================================
#Dialog for editing mesh settings
#===================================================================================================================================

class meshSettingsDialog(pyfrp_gui_basics.basicSettingsDialog):
	
	def __init__(self,mesh,parent):
		
		super(meshSettingsDialog,self).__init__(parent)
		
		self.mesh = mesh
		self.parent=parent
		self.nCharDisplayed=50
		
		#Labels
		self.lblVolSizePx = QtGui.QLabel("Element Size (px):", self)
		self.lblFromFile = QtGui.QLabel("From File?:", self)
		self.lblFnMesh = QtGui.QLabel("Filename mesh:", self)
		
		self.lblFnMeshVal = QtGui.QLabel("", self)
		self.updateFnMeshLbl()
		
		#LineEdits
		self.qleVolSizePx = QtGui.QLineEdit(str(self.mesh.volSizePx))
		
		self.doubleValid=QtGui.QDoubleValidator()
		self.intValid=QtGui.QIntValidator()
		
		self.qleVolSizePx.setValidator(self.doubleValid)
		
		self.qleVolSizePx.editingFinished.connect(self.setVolSizePx)
		
		#CheckBox
		self.cbFromFile = QtGui.QCheckBox('', self)
		self.updateCBs()
		self.connect(self.cbFromFile, QtCore.SIGNAL('stateChanged(int)'), self.checkFromFile)
		
		#Buttons
		self.btnFnMesh=QtGui.QPushButton('Change')
		
		self.btnFnMesh.connect(self.btnFnMesh, QtCore.SIGNAL('clicked()'), self.setFnMesh)
		
		#Layout
		self.fnMeshGrid = QtGui.QGridLayout()
		self.fnMeshGrid.addWidget(self.lblFnMeshVal,1,1)
		self.fnMeshGrid.addWidget(self.btnFnMesh,1,2)
		self.fnMeshGrid.setColumnMinimumWidth(1,150)
		
		self.grid.addWidget(self.lblVolSizePx,1,1)
		self.grid.addWidget(self.lblFromFile,2,1)
		self.grid.addWidget(self.lblFnMesh,3,1)
		
		self.grid.addWidget(self.qleVolSizePx,1,2)
		self.grid.addWidget(self.cbFromFile,2,2)
		self.grid.addLayout(self.fnMeshGrid,3,2)
		
		self.setWindowTitle("Mesh Settings")
		
		self.show()
		
	def setVolSizePx(self):
		self.mesh.setVolSizePx(float(str(self.qleVolSizePx.text())))
		
	def setFnMesh(self):
		
		fn = str(QtGui.QFileDialog.getExistingFile(self, "Select Mesh File",  self.parent.lastopen,"*.mol",))
		if fn=='':
			return
		
		self.mesh.setFnMesh(fn)
		
		self.parent.lastopen=folder
		self.updateFnMeshLbl()
		
	def updateFnMeshLbl(self):
		self.lblFnMeshVal.setText("..."+self.mesh.fnMesh[-self.nCharDisplayed:])	
		
	def checkFromFile(self,val):
		self.mesh.setFromFile(bool(2*val))

	def updateCBs(self):
		self.cbFromFile.setCheckState(2*int(self.mesh.fromFile))
	
#===================================================================================================================================
#Dialogs for mesh generation progress
#===================================================================================================================================

class genMeshProgressDialog(pyfrp_gui_basics.waitDialog):
	
	def __init__(self,parent):
		super(genMeshProgressDialog,self).__init__(parent)
		
		#Labels
		self.lblName.setText("Mesh generation in progress...")
		
		#Window title
		self.setWindowTitle('Mesh generation')
		    
		self.show()	

class genMeshThread(pyfrp_gui_basics.pyfrpThread):
	
	def __init__(self, embryo=None, parent=None):
		
		super(genMeshThread,self).__init__(parent)
		self.obj=embryo
		self.embryo=embryo
				
	def runTask(self,debug=False):
		self.embryo.simulation.mesh.genMesh()
		
#===================================================================================================================================
#Dialogs for mesh refine progress
#===================================================================================================================================

class refineMeshProgressDialog(pyfrp_gui_basics.waitDialog):
	
	def __init__(self,parent):
		super(refineMeshProgressDialog,self).__init__(parent)
		
		#Labels
		self.lblName.setText("Mesh refinement in progress...")
		
		#Window title
		self.setWindowTitle('Mesh refinement')
		    
		self.show()	

class refineMeshThread(pyfrp_gui_basics.pyfrpThread):
	
	def __init__(self, embryo=None, parent=None):
		
		super(refineMeshThread,self).__init__(parent)
		self.obj=embryo
		self.embryo=embryo
				
	def runTask(self,debug=False):
		self.embryo.simulation.mesh.refine()


		
#===================================================================================================================================
#Dialogs for forcing mesh density
#===================================================================================================================================

class basicForceMeshSettingsDialog(pyfrp_gui_basics.basicSettingsDialog):
	
	def __init__(self,mesh,parent):
		
		super(basicForceMeshSettingsDialog,self).__init__(parent)
	
		#Pass mesh
		self.mesh = mesh
		self.parent=parent
		self.embryo=self.mesh.simulation.embryo
		
		#Set default variables
		self.roiUsed=None
		self.debug=False
		self.findIdxs=True
		
		#Labels
		self.lblROI = QtGui.QLabel("ROI used:", self)
		self.lblDebug = QtGui.QLabel("Print debugging output:", self)
		self.lblFindIdxs = QtGui.QLabel("Update Indices?:", self)
		
		#LineEdits
		self.doubleValid=QtGui.QDoubleValidator()
		self.intValid=QtGui.QIntValidator()
		
		#Combobox
		self.comboROI = QtGui.QComboBox(self)
		self.updateROICombo()
		
		self.comboROI.activated[str].connect(self.setROI)   
		
		#CheckBox
		self.cbDebug = QtGui.QCheckBox('', self)
		self.cbFindIdxs = QtGui.QCheckBox('', self)
		
		self.connect(self.cbDebug, QtCore.SIGNAL('stateChanged(int)'), self.checkDebug)
		self.connect(self.cbFindIdxs, QtCore.SIGNAL('stateChanged(int)'), self.checkFindIdxs)
		
	
		#Layout	
		self.grid.addWidget(self.lblROI,1,1)
		self.grid.addWidget(self.lblFindIdxs,4,1)
		self.grid.addWidget(self.lblDebug,6,1)
		
		self.grid.addWidget(self.comboROI,1,2)
		self.grid.addWidget(self.cbFindIdxs,4,2)
		self.grid.addWidget(self.cbDebug,6,2)
		
		self.setWindowTitle("Basic Force Mesh Settings")
		
		self.show()
		
	def updateROICombo(self):
		for r in self.embryo.ROIs:
			self.comboROI.addItem(r.name)
		
		self.setROI()
		
	def setROI(self):
		self.roiUsed=self.mesh.simulation.embryo.ROIs[int(self.comboROI.currentIndex())]
	
	def checkDebug(self,val):
		self.debug=bool(2*val)
	
	def checkFindIdxs(self,val):
		self.findIdxs=bool(2*val)
		
	def donePressed(self):
		self.done(1)

class forceMeshSettingsDialog(basicForceMeshSettingsDialog):
	
	def __init__(self,mesh,parent):
		
		super(forceMeshSettingsDialog,self).__init__(mesh,parent)
	
		#Set default variables
		self.method='refine'
		self.stepPercentage=0.1
		self.maxCells=100000
		
		#Labels
		self.lblDensity = QtGui.QLabel("Desired Mesh Density:", self)
		self.lblStepPercentage = QtGui.QLabel("Step Size (%):", self)
		self.lblMethod = QtGui.QLabel("Method:", self)
		self.lblMaxCells = QtGui.QLabel("Maximum number of cells:", self)
		
		#LineEdits
		self.qleDensity = QtGui.QLineEdit("0")
		self.qleMaxCells = QtGui.QLineEdit(str(self.maxCells))
		self.qleStepPercentage = QtGui.QLineEdit(str(self.stepPercentage))
	
		self.qleDensity.setValidator(self.doubleValid)
		self.qleMaxCells.setValidator(self.intValid)
		self.qleStepPercentage.setValidator(self.doubleValid)
		
		self.qleDensity.editingFinished.connect(self.setDensity)
		self.qleMaxCells.editingFinished.connect(self.setMaxCells)
		self.qleStepPercentage.editingFinished.connect(self.setStepPercentage)
		
		#Combobox
		self.comboMethod = QtGui.QComboBox(self)
		self.addItem("refine")
		self.addWidget("method")
		
		self.comboMethod.activated[str].connect(self.setMethod)   
		
		#Layout	
		self.grid.addWidget(self.lblDensity,2,1)
		self.grid.addWidget(self.lblMethod,3,1)
		self.grid.addWidget(self.lblMaxCells,5,1)
		
		self.grid.addWidget(self.qleDensity,2,2)
		self.grid.addWidget(self.comboMethod,3,2)
		self.grid.addWidget(self.qleMaxCells,5,2)
	
		self.setWindowTitle("Force Global Mesh Density Settings")
		
		self.show()
	
	def updateROICombo(self):
		for r in self.embryo.ROIs:
			self.comboROI.addItem(r.name)
		
		self.setROI()
		self.qleDensity.setText(str(self.roiUsed.getMeshDensity()))
	
	def setDensity(self):
		self.density=float(str(self.qleDensity.text()))
	
	def setMaxCells(self):
		self.maxCells=int(str(self.qleMaxCells.text()))
		
	def setMethod(self,text):
		self.method=str(text)
	
	def setStepPercentage(self):
		self.stepPercentage=float(str(qleStepPercentage.text()))
	
	def getVals(self):
		return self.roiUsed,self.density,self.stepPercentage,self.debug,self.findIdxs,self.method,self.maxCells
	
	def donePressed(self):
		self.getVals()
		self.done(1)

#===================================================================================================================================
#Dialogs for forcing mesh density progress
#===================================================================================================================================

class refineROIMeshSettingsDialog(basicForceMeshSettingsDialog):
	
	def __init__(self,mesh,parent):
		
		super(refineROIMeshSettingsDialog,self).__init__(mesh,parent)
	
		#Set default variables
		self.factor=3.
		self.addZ=15
		
		#Labels
		self.lblFactor = QtGui.QLabel("Factor :", self)
		self.lblAddZ = QtGui.QLabel("Extend Box Mesh by z(px):", self)
		
		#LineEdits
		self.qleFactor = QtGui.QLineEdit(str(self.factor))
		self.qleAddZ = QtGui.QLineEdit(str(self.addZ))
	
		self.qleAddZ.setValidator(self.doubleValid)
		self.qleFactor.setValidator(self.doubleValid)
		
		self.qleFactor.editingFinished.connect(self.setFactor)
		self.qleAddZ.editingFinished.connect(self.setAddZ)
		
		#Layout
		self.grid.addWidget(self.lblFactor,2,1)
		self.grid.addWidget(self.lblAddZ,3,1)
		
		self.grid.addWidget(self.qleFactor,2,2)
		self.grid.addWidget(self.qleAddZ,3,2)
	
		self.setWindowTitle("Refine ROI Mesh Density Settings")
		
		self.show()
		
	def setAddZ(self):
		self.addZ=float(str(self.qleAddZ.text()))
		
	def setFactor(self):
		self.factor=float(str(self.qleFactor.text()))
	
	def getVals(self):
		return self.roiUsed,self.factor,self.addZ,self.debug,self.findIdxs
	
	def donePressed(self):
		if self.roiUsed==None:
			printWarning("No ROI selected. Going to do nothing.")
		else:
			self.roiUsed.refineInMeshByField(factor=self.factor,addZ=self.addZ,findIdxs=self.findIdxs,debug=self.debug)
		
		self.done(1)
	
			
#===================================================================================================================================
#Dialogs for forcing mesh density progress
#===================================================================================================================================

	
class forceMeshProgressDialog(pyfrp_gui_basics.waitDialog):
	
	def __init__(self,parent):
		super(forceMeshProgressDialog,self).__init__(parent)
		
		#Labels
		self.lblName.setText("Force mesh density in progress...")
		
		#Window title
		self.setWindowTitle('Mesh density enforcing')
		    
		self.show()	

class forceMeshThread(pyfrp_gui_basics.pyfrpThread):
	
	def __init__(self, embryo=None, parent=None):
		
		super(forceMeshThread,self).__init__(parent)
		self.obj=embryo
		self.embryo=embryo
		
	def runTask(self,roi,debug=False):
		self.embryo.simulation.mesh.forceMinMeshDensityInROI(roi)
		
				
		
			
			
		
		
				
		
		