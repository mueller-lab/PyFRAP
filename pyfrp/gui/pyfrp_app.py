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

#! /usr/bin/python

#=====================================================================================================================================
#Description
#=====================================================================================================================================

#PyFRAP is a free Python based software to analyze FRAP measurements. This file contains the main GUI helping to access the different 
#modules and classes:

#=====================================================================================================================================
#Importing Packages and Modules
#=====================================================================================================================================

#Standard packages
import sys
import os, os.path
import time
import copy as cpy
import functools
import code

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp_term import *
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_IO_module
from pyfrp.modules import pyfrp_plot_module

#PyFRAP GUIs
from pyfrp.gui import pyfrp_gui_molecule_dialogs
from pyfrp.gui import pyfrp_gui_embryo_dialogs
from pyfrp.gui import pyfrp_gui_ROI_manager
from pyfrp.gui import pyfrp_gui_geometry_dialogs
from pyfrp.gui import pyfrp_gui_gmsh_editor
from pyfrp.gui import pyfrp_gui_analysis_dialogs
from pyfrp.gui import pyfrp_gui_simulation_dialogs
from pyfrp.gui import pyfrp_gui_mesh_dialogs
from pyfrp.gui import pyfrp_gui_fit_dialogs
from pyfrp.gui import pyfrp_gui_pinning_dialogs
from pyfrp.gui import pyfrp_gui_statistics_dialogs
from pyfrp.gui import pyfrp_gui_basics

#PyFRAP Classes
from pyfrp.subclasses import pyfrp_conf
from pyfrp.subclasses import pyfrp_molecule

#QT
from PyQt4 import QtGui, QtCore

#Matplotlib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

"""
Apparently the NavigationToolbar naming has changed in newer matplotlib versions, thus
we need to test out some cases.
"""

try:
	from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as  NavigationToolbar
except ImportError:
	try:
		from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as  NavigationToolbar
	except ImportError:
		printWarning("Cannot import NavigationToolbar.")
		
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

#VTK 
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

#=====================================================================================================================================
#Main Simulation window
#=====================================================================================================================================

class pyfrp(QtGui.QMainWindow):
	
	#Initializes main GUI window
	def __init__(self, parent=None,redirect=False):
		QtGui.QWidget.__init__(self, parent)
		
		#-------------------------------------------
		#Version
		#-------------------------------------------
		
		from pyfrp import __version__ as version
		
		#-------------------------------------------
		#Initializing window
		#-------------------------------------------
		
		self.setWindowTitle('PyFRAP')
		self.setMinimumSize(400,300) 
		self.resize(1500,1000)
		self.dpi = 100
		self.version=version
		self.website="http://www.fml.tuebingen.mpg.de/de/mueller-group/software.html"
		self.pyfrpDir=pyfrp_misc_module.getConfDir().replace("configurations/","")
		
		#-------------------------------------------
		#Statusbar
		#-------------------------------------------
		
		self.statusBar().showMessage("Idle")
		
		#-------------------------------------------
		#Some variables
		#-------------------------------------------
		
		#Keep track of plotting tabs
		self.tabAxes=[]
		self.tabFigs=[]
		
		#Keep track of all molecules currently open
		self.molecules=[]
		
		#Keep track of which objects are currently used
		self.currMolecule=None
		self.currNode=None
		self.currMoleculeNode=None
		self.currObj=None
		
		#Keep track of which folder was recently used
		self.lastopen=os.getcwd()
	
		#-------------------------------------------
		#Menubar
		#-------------------------------------------
		
		self.menubar = self.menuBar()
		
		#Main Menubar Entries
		self.initMainMenubar()
		
		#File Menubar Entries
		self.initFileMenubar()
		
		#Edit Menubar Entries
		self.initEditMenubar()
		
		#Embryo Menubar Entries
		self.initEmbryoMenubar()
		
		#Analysis Menubar Entries
		self.initAnalysisMenubar()
		
		#Simulation Menubar Entries
		self.initSimulationMenubar()
		
		#Pinning Menubar Entries
		self.initPinningMenubar()
		
		#Fitting Menubar Entries
		self.initFittingMenubar()
		
		#Stats Menubar Entries
		self.initStatsMenubar()
	
		#-------------------------------------------
		#Embryo list
		#-------------------------------------------

		self.objectBar=QtGui.QTreeWidget()
		self.objectBar.setHeaderLabels(["Object","Analyzed","Simulated","Fitted"])
		self.objectBar.setColumnWidth(0,200)
		self.objectBar.setColumnWidth(1,75)
		self.objectBar.setColumnWidth(2,75)
		self.objectBar.itemClicked.connect(self.updatePropBar) 
		
		#-------------------------------------------
		#Property bar
		#-------------------------------------------
		
		self.propBar=QtGui.QListWidget()
		#self.propBar.itemDoubleClicked.connect(self.edit_prop)
		
		#-------------------------------------------
		#Console
		#-------------------------------------------
			
		self.console = PyInterp(self,redirect=redirect)
		self.console.initInterpreter(locals())
		
		#-------------------------------------------
		#Splitter
		#-------------------------------------------
		
		#Creating splitters
		self.horizontalSplitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
		self.verticalSplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
		
		#-------------------------------------------
		#Frames
		#-------------------------------------------
		
		#Creating frames for widgets
		self.objectBarFrame   = QtGui.QFrame()
		self.objectBarFrame   = self.createEmptyFrame(self.objectBarFrame)
		
		self.propBarFrame   = QtGui.QFrame()
		self.propBarFrame   = self.createEmptyFrame(self.propBarFrame)
		
		self.terminalFrame   = QtGui.QFrame()
		self.terminalFrame   = self.createEmptyFrame(self.terminalFrame)
		
		#-------------------------------------------
		#Plotting tabs
		#-------------------------------------------
		
		self.plotTabs = QtGui.QTabWidget()
		self.plotTabs.setTabsClosable(False)
		self.plotTabs.currentChanged.connect(self.currentTabChanged)
		self.plotTabs.tabCloseRequested.connect(self.currentTabClosed)
		
		#-------------------------------------------
		#Final Layout
		#-------------------------------------------
		
		##Adding widgets to splitters	
		self.horizontalSplitter.addWidget(self.objectBar)
		self.horizontalSplitter.addWidget(self.plotTabs)
		self.horizontalSplitter.addWidget(self.propBar)
		self.verticalSplitter.addWidget(self.horizontalSplitter)
		self.verticalSplitter.addWidget(self.console)
		
		#Setting default sizes of splitters
		self.verticalSplitter.setSizes([750,250])
		self.horizontalSplitter.setSizes([350,900,250])
		
		#Connecting splitter movement to figure size adjustment
		self.horizontalSplitter.splitterMoved.connect(self.adjustCanvas)
		self.verticalSplitter.splitterMoved.connect(self.adjustCanvas)
		
		#Create first plot tab
		self.createDummpyTab()
		
		#Load config file
		self.initConfiguration()
		
		self.setCentralWidget(self.verticalSplitter)
		QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Cleanlooks'))
		
		self.show()
		
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#List of Methods
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Initialization Methods
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Initiates main menubar
	
	def initMainMenubar(self):
		
		self.mbFile = self.menubar.addMenu('&File')
		self.mbEdit = self.menubar.addMenu('&Edit')
		self.mbEmbryo = self.menubar.addMenu('&Embryo')
		self.mbAnalysis = self.menubar.addMenu('&Data Analysis')
		self.mbSimulation = self.menubar.addMenu('&Simulation')
		self.mbPinning = self.menubar.addMenu('&Pinning')
		self.mbFitting = self.menubar.addMenu('&Fitting')
		self.mbStatistics = self.menubar.addMenu('&Statistics')
		
		return 
	def initFileMenubar(self):
		
		"""Creates entries of file menubar and connects actions with gui methods.
		"""
		
		newMoleculeButton = QtGui.QAction('New Molecule', self)
		self.connect(newMoleculeButton, QtCore.SIGNAL('triggered()'), self.newMolecule)
		
		loadMoleculeButton = QtGui.QAction('Open Molecule', self)
		self.connect(loadMoleculeButton, QtCore.SIGNAL('triggered()'), self.loadMolecule)
		
		saveMoleculeButton = QtGui.QAction('Save Molecule', self)
		self.connect(saveMoleculeButton, QtCore.SIGNAL('triggered()'), self.saveMolecule)
		
		exitButton = QtGui.QAction('Exit', self)
		exitButton.setShortcut('Ctrl+Q')	
		self.connect(exitButton, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))
	
	
		self.mbFile.addAction(newMoleculeButton)
		self.mbFile.addAction(loadMoleculeButton)
		self.recentMB=self.mbFile.addMenu('&Open recent')
		self.mbFile.addAction(saveMoleculeButton)
		
		self.mbFile.addAction(exitButton)
		
		return
	
	def initEditMenubar(self):
		
		"""Creates entries of edit menubar and connects actions with gui methods.
		"""
		
		editMoleculeButton = QtGui.QAction('Edit Molecule', self)
		self.connect(editMoleculeButton, QtCore.SIGNAL('triggered()'), self.editMolecule)
		
		removeMoleculeButton = QtGui.QAction('Remove Molecule', self)
		self.connect(removeMoleculeButton, QtCore.SIGNAL('triggered()'), self.removeMolecule)
		
		openWizardButton = QtGui.QAction('PyFRAP Wizard', self)
		self.connect(openWizardButton, QtCore.SIGNAL('triggered()'), self.openWizard)
		
		self.mbEdit.addAction(editMoleculeButton)
		self.mbEdit.addAction(removeMoleculeButton)
		self.mbEdit.addAction(openWizardButton)
		
	def initEmbryoMenubar(self):
		
		"""Creates entries of embryo menubar and connects actions with gui methods.
		"""
		
		newEmbryoButton = QtGui.QAction('New Embryo', self)
		self.connect(newEmbryoButton, QtCore.SIGNAL('triggered()'), self.newEmbryo)
		
		removeEmbryoButton = QtGui.QAction('Remove Embryo', self)
		self.connect(removeEmbryoButton, QtCore.SIGNAL('triggered()'), self.removeEmbryo)
		
		editEmbryoButton = QtGui.QAction('Edit Embryo', self)
		self.connect(editEmbryoButton, QtCore.SIGNAL('triggered()'), self.editEmbryo)
		
		loadEmbryoButton = QtGui.QAction('Load Embryo', self)
		self.connect(loadEmbryoButton, QtCore.SIGNAL('triggered()'), self.loadEmbryo)
		
		selectGeometryButton = QtGui.QAction('Select Geometry', self)
		self.connect(selectGeometryButton, QtCore.SIGNAL('triggered()'), self.selectGeometry)
		
		editGeometryButton = QtGui.QAction('Edit Geometry', self)
		self.connect(editGeometryButton, QtCore.SIGNAL('triggered()'), self.editGeometry)
		
		editGeoFileButton = QtGui.QAction('Edit GeoFile', self)
		self.connect(editGeoFileButton, QtCore.SIGNAL('triggered()'), self.editGeoFile)
		
		drawGeometryButton = QtGui.QAction('Plot Geometry', self)
		self.connect(drawGeometryButton, QtCore.SIGNAL('triggered()'), self.drawGeometry)
		
		roiManagerButton = QtGui.QAction('ROI Manager', self)
		self.connect(roiManagerButton, QtCore.SIGNAL('triggered()'), self.openROIManager)
		
		defaultROIButton = QtGui.QAction('Create Default ROIs', self)
		self.connect(defaultROIButton, QtCore.SIGNAL('triggered()'), self.createDefaultROIs)
		
		defaultROIWizardButton = QtGui.QAction('Default ROIs Wizard', self)
		self.connect(defaultROIWizardButton, QtCore.SIGNAL('triggered()'), self.defaultROIsWizard)
		
		indexROIButton = QtGui.QAction('Update ROIs indices', self)
		self.connect(indexROIButton, QtCore.SIGNAL('triggered()'), self.updateROIIdxs)
		
		self.mbEmbryo.addAction(newEmbryoButton)
		self.mbEmbryo.addAction(editEmbryoButton)
		self.mbEmbryo.addAction(loadEmbryoButton)
		self.mbEmbryo.addAction(removeEmbryoButton)
		
		self.geometryMB=self.mbEmbryo.addMenu('&Geometry')
		self.geometryMB.addAction(selectGeometryButton)
		self.geometryMB.addAction(editGeometryButton)
		self.geometryMB.addAction(editGeoFileButton)
		self.geometryMB.addAction(drawGeometryButton)
		
		self.roiMB=self.mbEmbryo.addMenu('&ROIs')
		self.roiMB.addAction(roiManagerButton)
		self.roiMB.addAction(defaultROIButton)
		self.roiMB.addAction(defaultROIWizardButton)
		
		self.roiMB.addAction(indexROIButton)
		
		return 
		
	def initAnalysisMenubar(self):
		
		"""Creates entries of analysis menubar and connects actions with gui methods.
		"""
		
		editAnalysisButton = QtGui.QAction('Analysis Settings', self)
		self.connect(editAnalysisButton, QtCore.SIGNAL('triggered()'), self.editAnalysis)
		
		analyzeEmbryoButton = QtGui.QAction('Analyze Embryo', self)
		self.connect(analyzeEmbryoButton, QtCore.SIGNAL('triggered()'), self.analyzeEmbryo)
		
		plotEmbryoAnalysisButton = QtGui.QAction('Plot Analysis Result of Embryo', self)
		self.connect(plotEmbryoAnalysisButton, QtCore.SIGNAL('triggered()'), self.plotAllDataTSOfEmbryo)
		
		self.mbAnalysis.addAction(editAnalysisButton)
		self.runAnalysisMB=self.mbAnalysis.addMenu('&Run Analysis')
		self.runAnalysisMB.addAction(analyzeEmbryoButton)
		
		
		self.plotAnalysisMB=self.mbAnalysis.addMenu('&Plotting')
		self.plotAnalysisMB.addAction(plotEmbryoAnalysisButton)
		
		return
	
	def initSimulationMenubar(self):
		
		"""Creates entries of simulation menubar and connects actions with gui methods.
		"""
		
		editSimulationButton = QtGui.QAction('Simulation Settings', self)
		self.connect(editSimulationButton, QtCore.SIGNAL('triggered()'), self.editSimulation)
		
		simulateEmbryoButton = QtGui.QAction('Simulate Embryo', self)
		self.connect(simulateEmbryoButton, QtCore.SIGNAL('triggered()'), self.simulateEmbryo)
		
		editMeshButton = QtGui.QAction('Mesh Settings', self)
		self.connect(editMeshButton, QtCore.SIGNAL('triggered()'), self.editMesh)
		
		genMeshButton = QtGui.QAction('Generate Mesh', self)
		self.connect(genMeshButton, QtCore.SIGNAL('triggered()'), self.generateMesh)
		
		refineMeshButton = QtGui.QAction('Refine Mesh', self)
		self.connect(refineMeshButton, QtCore.SIGNAL('triggered()'), self.refineMesh)
		
		
		forceMeshButton = QtGui.QAction('Force Mesh Density', self)
		self.connect(forceMeshButton, QtCore.SIGNAL('triggered()'), self.forceMeshDensity)
		
		forceROIMeshButton = QtGui.QAction('Refine Mesh in ROI', self)
		self.connect(forceROIMeshButton, QtCore.SIGNAL('triggered()'), self.forceROIMeshDensity)
		
		printMeshButton = QtGui.QAction('Mesh Stats', self)
		self.connect(printMeshButton, QtCore.SIGNAL('triggered()'), self.printMeshStats)
		
		plotMeshButton = QtGui.QAction('Plot Mesh', self)
		self.connect(plotMeshButton, QtCore.SIGNAL('triggered()'), self.plotMesh)
		
		plotMeshDensityButton = QtGui.QAction('Plot Mesh Density', self)
		self.connect(plotMeshDensityButton, QtCore.SIGNAL('triggered()'), self.plotMeshDensity)
		
		plotSimButton = QtGui.QAction('Plot Simulation', self)
		self.connect(plotSimButton, QtCore.SIGNAL('triggered()'), self.plotAllSimTSOfEmbryo)
		
		plotSimAndDataButton = QtGui.QAction('Plot Simulaten & Data', self)
		self.connect(plotSimAndDataButton, QtCore.SIGNAL('triggered()'), self.plotAllSimAndDataTSOfEmbryo)
		
		
		self.mbSimulation.addAction(editSimulationButton)
		self.runSimulationMB=self.mbSimulation.addMenu('&Run Simulation')
		self.runSimulationMB.addAction(simulateEmbryoButton)
		
		self.meshSimulationMB=self.mbSimulation.addMenu('&Mesh')
		self.meshSimulationMB.addAction(editMeshButton)
		self.meshSimulationMB.addAction(genMeshButton)
		self.meshSimulationMB.addAction(refineMeshButton)
		self.meshSimulationMB.addAction(forceMeshButton)
		self.meshSimulationMB.addAction(forceROIMeshButton)
		self.meshSimulationMB.addAction(printMeshButton)
		self.meshSimulationMB.addAction(plotMeshButton)
		self.meshSimulationMB.addAction(plotMeshDensityButton)
		
		
		self.plotSimulationMB=self.mbSimulation.addMenu('&Plotting')
		self.plotSimulationMB.addAction(plotSimButton)
		self.plotSimulationMB.addAction(plotSimAndDataButton)

	def initPinningMenubar(self):
		
		"""Creates entries of pinning menubar and connects actions with gui methods.
		"""
		
		defaultPinButton = QtGui.QAction('Default Pinning', self)
		self.connect(defaultPinButton, QtCore.SIGNAL('triggered()'), self.defaultPinEmbryo)
		
		idealPinButton = QtGui.QAction('Ideal Pinning', self)
		self.connect(idealPinButton, QtCore.SIGNAL('triggered()'), self.idealPinEmbryo)
		
		self.mbPinning.addAction(defaultPinButton)
		self.mbPinning.addAction(idealPinButton)
		
	def initFittingMenubar(self):
		
		"""Creates entries of fitting menubar and connects actions with gui methods.
		"""
		
		newFitButton = QtGui.QAction('New fit', self)
		self.connect(newFitButton, QtCore.SIGNAL('triggered()'), self.newFit)
		
		removeFitButton = QtGui.QAction('Remove fit', self)
		self.connect(removeFitButton, QtCore.SIGNAL('triggered()'), self.removeFit)
		
		editFitButton = QtGui.QAction('Edit fit', self)
		self.connect(editFitButton, QtCore.SIGNAL('triggered()'), self.editFit)
		
		#editmultfit = QtGui.QAction('Edit multiple fits', self)
		#self.connect(editmultfit, QtCore.SIGNAL('triggered()'), self.edit_mult_fit)
		
		#copyfit = QtGui.QAction('Copy fit', self)
		#self.connect(copyfit, QtCore.SIGNAL('triggered()'), self.copy_fit)
		
		#copyfitforall = QtGui.QAction('Copy fit into all embryos', self)
		#self.connect(copyfitforall, QtCore.SIGNAL('triggered()'), self.copy_fit_to_all)
		
		performFitButton = QtGui.QAction('Perform fit', self)
		self.connect(performFitButton, QtCore.SIGNAL('triggered()'), self.performFit)
		
		
		#fitall = QtGui.QAction('Perform all fits in molecule', self)
		#self.connect(fitall, QtCore.SIGNAL('triggered()'), self.perform_fits_molecule)
		
		plotFitButton = QtGui.QAction('Plot fit', self)
		self.connect(plotFitButton, QtCore.SIGNAL('triggered()'), self.plotFit)
		
		#plottrackfit = QtGui.QAction('Plot fitting progress', self)
		#self.connect(plottrackfit, QtCore.SIGNAL('triggered()'), self.plot_track_fit)
		
		printFitButton = QtGui.QAction('Print fit results', self)
		self.connect(printFitButton, QtCore.SIGNAL('triggered()'), self.printFitResults)
		
		
		self.mbFitting.addAction(newFitButton)
		self.mbFitting.addAction(editFitButton)
		self.mbFitting.addAction(removeFitButton)
		self.mbFitting.addAction(performFitButton)
		self.mbFitting.addAction(printFitButton)
		
		
		self.plotFittingMB=self.mbFitting.addMenu('&Plotting')
		self.plotFittingMB.addAction(plotFitButton)
		
		return
		
	def initStatsMenubar(self):
		
		"""Creates entries of statistics menubar and connects actions with gui methods.
		"""
		
		selectFitsButton = QtGui.QAction('Select Fits', self)
		self.connect(selectFitsButton, QtCore.SIGNAL('triggered()'), self.selectFits)
		
		selectCrucialParametersButton = QtGui.QAction('Select CrucialParameters', self)
		self.connect(selectCrucialParametersButton, QtCore.SIGNAL('triggered()'), self.selectCrucialParameters)
		
		
		summarizeMoleculeButton = QtGui.QAction('Summarize Molecule', self)
		self.connect(summarizeMoleculeButton, QtCore.SIGNAL('triggered()'), self.summarizeMolecule)
		
		
		self.mbStatistics.addAction(selectFitsButton)
		self.mbStatistics.addAction(selectCrucialParametersButton)
		self.mbStatistics.addAction(summarizeMoleculeButton)
	
	def closeEvent(self, event):
		
		"""Closes GUI and saves configuration.
		
		Args:
			event (QtCore.closeEvent): Close event triggered by close signal.
		
		"""
		
		reply = QtGui.QMessageBox.question(self, 'Message',"Are you sure you want to quit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
	
		if reply == QtGui.QMessageBox.Yes:
			fn=pyfrp_misc_module.getConfDir()+"lastConfiguration.conf"
			self.config.consoleHistory=self.console.history
			self.config.save(fn=fn)
			event.accept()
		else:
			event.ignore()
		
		return
		
	#def show_about(self):
		
		#ret=pyfrp_subwin.about_dialog(self).exec_()
	
	
	def createEmptyFrame(self,frame):
		
		"""Creates frame around widget.
		
		Args:
			frame (QtGui.QWidget): Widget to be framed
		
		Returns: 
			QtGui.QWidget: Framed Widget
		
		"""
			
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		frame.setBackgroundRole(QtGui.QPalette.Light)
		frame.setAutoFillBackground(True)        
		frame.setLineWidth(1)
		frame.setFrameShadow(frame.Sunken)
		
		return frame
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Configuration handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Initialize Configuration from file
	
	def initConfiguration(self):
		
		fn=pyfrp_misc_module.getConfDir()+"lastConfiguration.conf"
		
		if os.path.isfile(fn):
			self.config=pyfrp_IO_module.loadFromPickle(fn)
		else:
			self.config=pyfrp_conf.configuration()
		
		self.console.history=list(self.config.consoleHistory)
		
		self.updateConfig()
		
		return self.config
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Update configuration of GUI to match with config
	
	def updateConfig(self):
		
		self.console.setHidden(self.config.termHidden)
		self.propBar.setHidden(self.config.propHidden)
		self.plotTabs.setHidden(self.config.plotHidden)
		
		self.verticalSplitter.refresh()
		self.horizontalSplitter.refresh()
		#self.adjustCanvas()
		
		self.updateRecentMBs()
		
		return self.config
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Makes filename the most recently opened file, adds it to config.recentFiles and updates menubar  
	
	def appendRecent(self,fn):
	
		self.config.addRecentFile(fn)
		self.recentMB.clear()
		self.updateRecentMBs()
		
		return self.recentMB
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Updates recently opened menubar  
	
	def updateRecentMBs(self):
		self.recentActions=[]
		
		for i in range(len(self.config.recentFiles)): 
			if i>5:
				self.config.recentFiles.pop(i)
			else:
				self.recentActions.append(QtGui.QAction(self.config.recentFiles[i], self))
				item=self.recentActions[i]
				item.triggered.connect(functools.partial(self.openMolecule,self.config.recentFiles[i]))
				
				self.recentMB.addAction(item)	
				
		return self.recentMB		
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Object Bar handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Highlists the currently used Node of the objectBar
	
	def highlightCurrentObjectNode(self):
		self.objectBar.setCurrentItem(self.currNode)
		return self.objectBar
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Checks if the current selected node is of certain type and returns respective error popups
	
	def checkSelectedNode(self,typ="any"):
		
		if self.objectBar.currentItem()==None:
			QtGui.QMessageBox.critical(None, "Error","Nothing selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return False
		if typ=="any":
			return True
		if typ==self.currNodeType:
			return True
		
		
		return False
				
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Adds Embryo to ObjectBar

	def embryo2ObjectBar(self,embryo,parentNode):
		
		#Embryo Status
		analyzed=str(int(embryo.isAnalyzed()))
		simulated=str(int(embryo.isSimulated()))
		fitted=str(int(embryo.isFitted()))
		
		#Create new node
		newEmbryoNode=QtGui.QTreeWidgetItem(parentNode,[embryo.name,analyzed,simulated,fitted])
		
		#Geometry node
		newGeometryNode=QtGui.QTreeWidgetItem(newEmbryoNode,['Geometry','','',''])
		
		#ROIs node
		newROIsNode=QtGui.QTreeWidgetItem(newEmbryoNode,['ROIs','','',''])
		
		#Add ROIs
		for roi in embryo.ROIs:
			newROINode=QtGui.QTreeWidgetItem(newROIsNode,[roi.name,str(int(roi.isAnalyzed())),str(int(roi.isSimulated())),str(int(roi.isFitted()))])
		
		#Analysis node
		newAnalysisNode=QtGui.QTreeWidgetItem(newEmbryoNode,['Analysis',analyzed,'',''])
		
		#Simulation node
		newSimulationNode=QtGui.QTreeWidgetItem(newEmbryoNode,['Simulation','',simulated,''])
		
		#if embryo.simulation!=None:
		newMeshNode=QtGui.QTreeWidgetItem(newSimulationNode,['Mesh','','',''])
			
		#Fits node
		newFitsNode=QtGui.QTreeWidgetItem(newEmbryoNode,['Fits','','',''])
		
		#Add Fits
		for fit in embryo.fits:
			newFitNode=QtGui.QTreeWidgetItem(newFitsNode,[fit.name,'','',str(int(fit.isFitted()))])
		
		#Expand Fits node
		self.objectBar.expandItem(newFitsNode)	
		
		return newEmbryoNode		
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Update ObjectBar Properties of embryo Node

	def updateEmbryoNodeProps(self,embryoNode):
		currEmbryoNode = self.getCurrentEmbryoNode()
		currEmbryo = self.getCurrentEmbryo()
		
		currEmbryoNode.setText(1,str(int(currEmbryo.isAnalyzed())))
		currEmbryoNode.setText(2,str(int(currEmbryo.isSimulated())))
		currEmbryoNode.setText(3,str(int(currEmbryo.isFitted())))
		
		return
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Update ObjectBar Properties of embryo Node

	def updateEmbryoNodeChildren(self):
		self.updateFitsNodeChildren()
		self.updateROIsNodeChildren()
		return
	
	def getChildByName(self,node,name):
		for i in range(node.childCount()):
			if name==str(node.child(i).data(0,0).toString()):
				return node.child(i)
		return None
	
	def removeAllNodeChildren(self,node):
		nChild=node.childCount()
		for i in range(nChild):
			node.takeChild(0)
			
		return node

	def updateFitsNodeChildren(self):
		
		currEmbryo=self.getCurrentEmbryo()
		currEmbryoNode=self.getCurrentEmbryoNode()
		if currEmbryoNode!=None:
			fitsNode=self.getChildByName(currEmbryoNode,'Fits')
			self.removeAllNodeChildren(fitsNode)
			for fit in currEmbryo.fits:
				newFitNode=QtGui.QTreeWidgetItem(fitsNode,[fit.name,'','',str(int(fit.isFitted()))])
		self.getCurrentObjNode()
		return newFitNode
		
		
	def updateROIsNodeChildren(self):
		
		currEmbryo=self.getCurrentEmbryo()
		currEmbryoNode=self.getCurrentEmbryoNode()
		if currEmbryoNode!=None:
			roisNode=self.getChildByName(currEmbryoNode,'ROIs')
			self.removeAllNodeChildren(roisNode)
			for roi in currEmbryo.ROIs:
				newROINode=QtGui.QTreeWidgetItem(roisNode,[roi.name,str(int(roi.isAnalyzed())),str(int(roi.isSimulated())),str(int(roi.isFitted()))])
		self.getCurrentObjNode()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Object Node from objectBar

	def getCurrentObjNode(self):
		
		self.currNode=self.objectBar.currentItem()
		
		if self.currNode.parent()==None:
			self.currNodeType='molecule'
		else:
			if self.currNode.parent().parent()==None:
				self.currNodeType='embryo'
			elif self.currNode.data(0,0).toString()=="Fits":	
				self.currNodeType='fits'
			elif self.currNode.parent().data(0,0).toString()=="Fits":	
				self.currNodeType='fit'
			elif self.currNode.data(0,0).toString()=="Simulation":	
				self.currNodeType='simulation'
			elif self.currNode.data(0,0).toString()=="Analysis":	
				self.currNodeType='analysis'
			elif self.currNode.data(0,0).toString()=="Mesh":	
				self.currNodeType='mesh'
			elif self.currNode.data(0,0).toString()=="Geometry":	
				self.currNodeType='geometry'
			elif self.currNode.data(0,0).toString()=="ROIs":	
				self.currNodeType='rois'
			elif self.currNode.parent().data(0,0).toString()=="ROIs":	
				self.currNodeType='roi'
			
		return self.currNode, self.currNodeType		
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Object selected in objectBar

	def getCurrentObj(self):
		
		#Find out current node and which type it is
		self.getCurrentObjNode()
		
		#Find current molecule node 
		self.getCurrentMoleculeNode()
		
		#Find corresponding molecule object
		self.getCurrentMolecule()
		
		if self.currNodeType=="molecule":
			self.currObj=self.currMolecule
		else:	
			currEmbryo=self.getCurrentEmbryo()
			
			if self.currNodeType=="embryo":
				self.currObj=currEmbryo
			elif self.currNodeType=="analysis":
				self.currObj=currEmbryo.analysis
			elif self.currNodeType=="geometry":
				self.currObj=currEmbryo.geometry
			elif self.currNodeType=="simulation":
				self.currObj=currEmbryo.simulation
			elif self.currNodeType=="mesh":
				self.currObj=currEmbryo.simulation.mesh
			elif self.currNodeType=="fit":
				self.currObj=currEmbryo.getFitByName(self.currNode.data(0,0).toString())
			elif self.currNodeType=="roi":
				self.currObj=currEmbryo.getROIByName(self.currNode.data(0,0).toString())
			else:
				self.currObj=None
				
		return self.currObj		
				
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Molecule Node
		
	def getCurrentMoleculeNode(self):
		if self.currNodeType=='molecule':
			self.currMoleculeNode=self.currNode
		elif self.currNodeType=='embryo':
			self.currMoleculeNode=self.currNode.parent()
		elif self.currNodeType in ["analysis","simulation","fits","geometry","rois"]:
			self.currMoleculeNode=self.currNode.parent().parent()
		elif self.currNodeType in ["mesh","roi","fit"]:
			self.currMoleculeNode=self.currNode.parent().parent().parent()
		
		
		###NOTE: Try later if this is the smarter way
		#ind=self.objectBar.indexOfTopLevelItem(self.currNode)
		#self.currMoleculeNode=self.objectBar.topLevelItem(ind)
		
		return self.currMoleculeNode
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Embryo Node
		
	def getCurrentEmbryoNode(self):
		if self.currNodeType=='molecule':
			currEmbryoNode=None
		elif self.currNodeType=='embryo':
			currEmbryoNode=self.currNode
		elif self.currNodeType in ["analysis","simulation","fits","geometry","rois"]:
			currEmbryoNode=self.currNode.parent()
		elif self.currNodeType in ["mesh","roi","fit"]:
			currEmbryoNode=self.currNode.parent().parent()
		
		return currEmbryoNode

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Molecule
	
	def getCurrentMolecule(self):	
		molNames=pyfrp_misc_module.objAttrToList(self.molecules,"name")
		self.currMolecule=self.molecules[molNames.index(self.currMoleculeNode.data(0,0).toString())]
		return self.currMolecule
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Embryo
	
	def getCurrentEmbryo(self):
		currEmbryoNode=self.getCurrentEmbryoNode()
		if currEmbryoNode==None:
			return None
		
		embryoNames=pyfrp_misc_module.objAttrToList(self.currMolecule.embryos,"name")
		currEmbryo=self.currMolecule.embryos[embryoNames.index(currEmbryoNode.data(0,0).toString())]
		return currEmbryo
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Porperty Bar handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Returns Current Molecule
	
	def showObjPropsInBar(self,obj,maxArraySize=3):	
		
		if obj==None:
			return self.propBar
		
		for item in vars(obj):
			
			###NOTE: Creates deprecation warning on Windows
			if isinstance(vars(obj)[str(item)],(int,float,str)) or vars(obj)[str(item)]==None:
				pass
			elif len(np.shape(vars(obj)[str(item)]))>0 and np.shape(vars(obj)[str(item)])[0]<maxArraySize:
				pass
			else:
				continue
			
			self.propBar.addItem(item+"="+str(vars(obj)[str(item)]))
		
		return self.propBar 
		
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Updates Property bar
		
	def updatePropBar(self):
		
		#Clearing property list
		self.propBar.clear()
		
		#Get Current Object
		self.getCurrentObj()
			
		#Show properties 
		self.showObjPropsInBar(self.currObj)
		
		#Sort prop_list
		self.propBar.sortItems()
		
		return self.propBar 

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Molecule handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Creates new molecule object and asks for name, then automatically adds it to objectBar
	
	def newMolecule(self):
		
		#Generate name
		moleculeNames=pyfrp_misc_module.objAttrToList(self.molecules,"name")
		newName=pyfrp_misc_module.enumeratedName("newMolecule",moleculeNames,sep="_")
		
		#Create new molecule object
		mol=pyfrp_molecule.molecule(newName)
		pyfrp_gui_molecule_dialogs.moleculeDialog(mol,self).exec_()
		
		#Make current molecule
		self.currMolecule=mol
		
		#Add new molecule to list of embryos
		self.molecules.append(mol)
		
		#Add to objectBar
		self.currNode=QtGui.QTreeWidgetItem(self.objectBar,[self.currMolecule.name])
		
		#Memorize molecule node
		self.currMoleculeNode=self.currNode
		
		#Highligth new node
		self.highlightCurrentObjectNode()
			
		#Show molecule properties
		self.updatePropBar()
		
		#Check if names are alright
		self.checkObjNames()
		
		return self.currMolecule
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Edit Molecule
	
	def editMolecule(self):
		
		self.getCurrentMolecule()
		ret = pyfrp_gui_molecule_dialogs.moleculeDialog(self.currMolecule,self).exec_()
		
		self.currMoleculeNode.setText(0,self.currMolecule.name)
		
		self.updatePropBar()
		
		return self.currMolecule
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Removes Molecule from PyFRAP
		
	def removeMolecule(self):
		
		#Check if molecule or subnode selected
		if not self.checkSelectedNode():
			return 
				
		#Remove from list of molecules used
		self.molecules.remove(self.currMolecule)
		
		#Remove from objectBar
		ind=self.objectBar.indexOfTopLevelItem(self.currMoleculeNode)
		self.objectBar.takeTopLevelItem(ind)
		
		#Update current Node
		self.currNode=self.objectBar.currentItem()
		
		#Update PropBar
		if self.currNode!=None:
			self.updatePropBar()
		else:
			self.propBar.clear()		
		
		return self.molecules
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Saves molecule object into pickle file
	
	def saveMolecule(self):
		
		#Check if molecule or subnode selected
		if not self.checkSelectedNode():
			return 
		
		#Filename dialog
		fn=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen+"/"+self.currMolecule.name+".mol","*.mol",)
		fn=str(fn)
			
		#Check if folder exists
		if not os.path.exists(os.path.dirname(fn)):
			return
		
		#Remember folder
		self.lastopen=os.path.dirname(str(fn))
			
		#Ask if extract save
		#reply = QtGui.QMessageBox.question(self, 'Message',"Do you want to save  (Will only keep essential results and delete everything else)?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
		
		#No slim save
		#if reply == QtGui.QMessageBox.No:
		self.currMolecule.save(fn=fn)
		#else:
			
			##Temporarily save molecule
			#if self.config.backup2File:
				#fn_backup=self.lastopen+"/"+self.currMolecule.name+"_backup.mol"
				#self.currMolecule.save(fn=fn_backup)
				#tempMolecule=self.currMolecule
			
			##Temporarily make a copy of molecule into memory
			#if self.config.backup2Memory:
				#tempMolecule=cpy.deepcopy(self.currMolecule)
			
			####NOTE: Need to specify what needs to be deleted for slim save
			##Make a function in embryo object for this
			##Delete big stuff
			##for temp_emb in temp_mol.embryos:
			
				###Deleting all image data to reduce file size
				##temp_emb.masks_embryo=[]
				
			#if self.config.backup2File:
				#self.currMolecule=pyfrp_IO_module.loadFromPickle(fn_backup)
				#os.remove(fn_backup)
			#if self.config.backup2Memory:
				#self.currMolecule=cpy.deepcopy(tempMolecule)
				#tempMolecule=None
				
		self.appendRecent(fn)
		
		return fn

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Load molecule
	
	def loadMolecule(self):
		
		fnLoad=QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.lastopen,"*.mol",)
		if fnLoad=='':
			return
		
		self.currMolecule=self.openMolecule(fnLoad)
		
		return self.currMolecule
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Open molecule
	
	def openMolecule(self,fnLoad):
		
		#Memorize last path and append to recently opened
		self.lastopen=os.path.dirname(str(fnLoad))
		self.appendRecent(fnLoad)
		
		#Load molecule object
		self.currMolecule=pyfrp_IO_module.loadFromPickle(fnLoad)
		
		#Update Version
		self.currMolecule.updateVersion()
		
		#Add molecule to list of molecules
		self.molecules.append(self.currMolecule)
		
		#Adding molecule to sidebar
		self.currNode=QtGui.QTreeWidgetItem(self.objectBar,[self.currMolecule.name,"","",""])
		self.currMoleculeNode=self.currNode
		
		for embryo in self.currMolecule.embryos:
			self.embryo2ObjectBar(embryo,self.currMoleculeNode)
						
		self.objectBar.expandItem(self.currNode)
		
		self.checkObjNames()
		
		return self.currMolecule
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Embryo handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	def newEmbryo(self):
		
		"""Creates new embryo object and adds it to molecule.
		
		Will first run :py:func:`selectEmbryoGenMode` to choose way of embryo generation.
		Will create embryo instance via :py:func:`LSM2EmbryoWizard` if selected, ohterwise
		will create default embryo.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: New embryo object.
		
		"""
		
		#Check if molecule or subnode selected
		if not self.checkSelectedNode():
			return 
		
		ret=self.selectEmbryoGenMode()
		if ret==-1:
			return
		elif ret==0:
			#Generate name
			embryoNames=pyfrp_misc_module.objAttrToList(self.currMolecule.embryos,"name")
			newName=pyfrp_misc_module.enumeratedName("newEmbryo",embryoNames,sep="_")
			
			#Create new embryo
			newEmbryo=self.currMolecule.newEmbryo(newName)
			
			#Append to object bar
			newNode=self.embryo2ObjectBar(newEmbryo,self.currMoleculeNode)
			self.objectBar.setCurrentItem(newNode)
			self.getCurrentObj()
			
		#Call embryo editor
		self.editEmbryo()
		
		#Geometry Editor
		self.selectGeometry()		
		self.editGeometry()
		
		return self.getCurrentEmbryo()
	
	def loadEmbryo(self):
		
		"""Loads saved embryo object and adds it to currently selected molecule file.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: Loaded embryo object.
		
		"""
		
		#Check if molecule or subnode selected
		if not self.checkSelectedNode():
			return 
		
		#Get filename
		fnLoad=QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.lastopen,"*.emb",)
		
		#Load and add to molecule
		newEmbryo=pyfrp_misc_module.loadFromPickle(fnLoad)
		self.currMolecule.addEmbryo(newEmbryo)
		
		#Append to object bar
		newNode=self.embryo2ObjectBar(newEmbryo,self.currMoleculeNode)
		self.objectBar.setCurrentItem(newNode)
		self.getCurrentObj()
		
		return newEmbryo
		
	def editEmbryo(self):
		
		"""Opens main embryo editing dialog.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: Edited embryo object.
		
		"""
		
		currEmbryoNode=self.getCurrentEmbryoNode()
		currEmbryo=self.getCurrentEmbryo()
		
		if self.checkSelectedNode(typ='molecule') or self.objectBar.currentItem()==None:
			return 
		
		ret=pyfrp_gui_embryo_dialogs.embryoDialog(currEmbryo,self).exec_()
		
		currEmbryoNode.setText(0,currEmbryo.name)
		
		self.updatePropBar()
			
		return currEmbryo
				
	def removeEmbryo(self):
		
		"""Removes selected embryo object both from molecule object as well as 
		main GUI.
		
		Returns:
			list: Updated pyfrp.subclasses.pyfrp_molecule.embryos list
		
		"""
		
		currEmbryoNode=self.getCurrentEmbryoNode()
		currEmbryo=self.getCurrentEmbryo()
		
		if self.checkSelectedNode(typ='molecule') or self.objectBar.currentItem()==None:
			return 
		
		#Remove from list of embryo objects
		self.currMolecule.embryos.remove(currEmbryo)
		
		#Remove from sidebar
		ind=self.currMoleculeNode.indexOfChild(currEmbryoNode)
		self.currMoleculeNode.takeChild(ind)
		
		self.objectBar.setCurrentItem(self.currMoleculeNode)
		
		self.updatePropBar()
		
		return self.currMolecule.embryos

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Checks if molecules or embryos have the same name

	def checkObjNames(self):
		
		b=False
		
		b=self.checkMolNames()
		
		for mol in self.molecules:
			b=self.checkEmbryoNames(mol=mol)
		
		return b
				
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Checks if molecules have the same name
			
	def checkMolNames(self):
		
		b=False
		
		molNames=pyfrp_misc_module.objAttrToList(self.molecules,"name")
		for name in molNames:
			if molNames.count(name)>1:
				b=True
				printWarning("Molecule with name " + name + " exists twice. This can lead to problems")
		
		return b
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Checks if embryos have the same name
			
	def checkEmbryoNames(self,mol=None):
		
		mol = pyfrp_misc_module.assignIfVal(mol,self.currMolecule,None)
		
		b=False
		
		embryoNames=pyfrp_misc_module.objAttrToList(self.currMolecule.embryos,"name")
		for name in embryoNames:
			if embryoNames.count(name)>1:
				b=True
				printWarning("Embryo with name " + name + " exists twice in molecule " + mol.name +". This can lead to problems")
		
		return b
	
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Wizards
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	def LSM2EmbryoWizard(self):
		
		"""Opens :py:class:`pyfrp.gui.pyfrp_gui_embryo_dialogs.lsmWizard` which lets
		one to extract and sort microscope files and then creates new 
		:py:class:`pyfrp.subclasses.pyfrp_embryo.embryo` including new data.
		
		.. note:: Uses Fiji. Fiji path needs to be set properly for success.
		
		"""
		
		lsmWizard = pyfrp_gui_embryo_dialogs.lsmWizard(self)
		if lsmWizard.exec_():
			embryo = lsmWizard.getEmbryo()
			if embryo==None:
				return	
			
			self.currMolecule.addEmbryo(embryo)
			
			newNode=self.embryo2ObjectBar(embryo,self.currMoleculeNode)
			self.objectBar.setCurrentItem(newNode)
			self.getCurrentObj()
					
	def openWizard(self):
		
		"""Wizard for PyFRAP analysis.
		
		Leads user through following steps:
		
			* Embryo creation :py:func:`newEmbryo`.
			* ROI creation :py:func:`createDefaultROIs`.
			* ROI handling :py:func:`openROIManager`.
			* Mesh generation :py:func:`generateMesh`.
			* Data analysis :py:func:`analyzeEmbryo`.
			* FRAP simulation :py:func:`simulateEmbryo`.
			* Fit creation :py:func:`newFit`.
			* Fit plotting :py:func:`plotFit`.
			
		"""	
		
		self.newEmbryo()
		self.defaultROIsWizard()
		self.openROIManager()
		self.editGeometry()
		self.getCurrentEmbryo().newSimulation()
		self.generateMesh()
		self.updateROIIdxs()
		self.analyzeEmbryo()
		self.simulateEmbryo()
		self.newFit()
		self.plotFit()
	
	def selectEmbryoGenMode(self):
		
		"""Opens :py:class:`pyfrp.gui.pyfrp_gui_embryo_dialogs.wizardSelector`
		and lets user select embryo generation mode. 
		
		Returns:
			int: 1 if user chose from microscope, 0 if not, -1 else.
		
		"""
		
		selectWizard = pyfrp_gui_embryo_dialogs.wizardSelector(self)
		if selectWizard.exec_():
			useLSM = selectWizard.getUseLSM()
			if useLSM==None:
				return -1
			
			if useLSM:
				self.LSM2EmbryoWizard()
				return 1
			return 0
				
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Embryo Menubar Methods
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	def openROIManager(self):
		
		"""Open ROI Manager for current embryo.
		"""
		
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			ret=pyfrp_gui_ROI_manager.ROImanager(currEmbryo,self).exec_()
			self.updateROIsNodeChildren()
		return
	
	def createDefaultROIs(self):
		
		"""Create default ROIs for current embryo.
		"""
		
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			ret=pyfrp_gui_ROI_manager.defaultROIsDialog(currEmbryo,self).exec_()
			self.updateROIsNodeChildren()
		return
	
	def defaultROIsWizard(self):
		
		"""Small Wizard that helps creating a set of default ROIs for current embryo.
		"""
		
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		#Select Slice ROI
		sliceSelector = pyfrp_gui_ROI_manager.ROISelector(currEmbryo,'radialSlice',self,msg='Select Slice ROI',sliceHeight=currEmbryo.sliceHeightPx,sliceWidth=currEmbryo.sliceWidthPx,color='g',name='Slice')
		if sliceSelector.exec_():
			sliceROI=sliceSelector.getROI()
			
		#Select Bleached ROI
		squSelector = pyfrp_gui_ROI_manager.ROISelector(currEmbryo,'squareSlice',self,msg='Select Bleached ROI',sliceHeight=currEmbryo.sliceHeightPx,sliceWidth=currEmbryo.sliceWidthPx,color='b',name='Bleached Square')
		if squSelector.exec_():
			bleachedROI=squSelector.getROI()
		
		#Select Rim ROI
		if hasattr(sliceROI,'radius'):
			radius=sliceROI.radius*0.66
			center=sliceROI.center
			rimSelector = pyfrp_gui_ROI_manager.ROISelector(currEmbryo,'radialSlice',self,msg='Select Rim ROI',radius=radius,center=center,sliceHeight=currEmbryo.sliceHeightPx,sliceWidth=currEmbryo.sliceWidthPx,color='y',name='Slice rim')
		else:
			radius=133.
			rimSelector = pyfrp_gui_ROI_manager.ROISelector(currEmbryo,'radialSlice',self,msg='Select Rim ROI',radius=radius,sliceHeight=currEmbryo.sliceHeightPx,sliceWidth=currEmbryo.sliceWidthPx,color='y',name='Slice rim')
			
		if rimSelector.exec_():
			rimROI=rimSelector.getROI()
		
		#Build default ROIs
		radius=300.
		center=[256.,256.]
		if hasattr(sliceROI,'radius'):
			radius=sliceROI.radius
			center=sliceROI.center
		else:
			if currEmbryo.geometry!=None:
				if hasattr(currEmbryo.geometry,'imagingRadius'):
					radius=currEmbryo.geometry.imagingRadius
				else:
					if hasattr(currEmbryo.geometry,'radius'):
						radius=currEmbryo.geometry.radius
		
		currEmbryo.genDefaultROIs(center,radius,rimFactor=0.66,masterROI=sliceROI,bleachedROI=bleachedROI,rimROI=rimROI,clean=True)
		
		self.updateROIsNodeChildren()

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Updates all ROI Idxs
		
	def updateROIIdxs(self):
		
		currEmbryo=self.getCurrentEmbryo()
	
		if currEmbryo!=None:
		
			#Genereate wait popup
			self.progressDialog=pyfrp_gui_ROI_manager.indexProgressDialog(None)
			
			#Make backup copy of embryo
			self.originalObj=currEmbryo
			self.backupObj=cpy.deepcopy(currEmbryo)
			
			self.statusBar().showMessage("Indexing ROIs of embryo  " + currEmbryo.name)
			
			#Generate worker and generate Qthread
			self.task=pyfrp_gui_basics.pyfrpThread()
			self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.computeROIIdxs,signal=self.task.progressSignal)
			
			#Start
			self.initTask()
			
		return

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Open Geometry Editor for current embryo
		
	def selectGeometry(self):
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo!=None:
			ret=pyfrp_gui_geometry_dialogs.geometrySelectDialog(currEmbryo,self).exec_()
		return	
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Open Geometry Editor for current embryo
		
	def editGeometry(self):
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo!=None:
			if "zebrafishDomeStage" in currEmbryo.geometry.typ:
				ret=pyfrp_gui_geometry_dialogs.zebrafishDomeStageDialog(currEmbryo.geometry,self).exec_()
			elif "cylinder" in currEmbryo.geometry.typ:
				ret=pyfrp_gui_geometry_dialogs.cylinderDialog(currEmbryo.geometry,self).exec_()
			elif "xenopusBall" in currEmbryo.geometry.typ:
				ret=pyfrp_gui_geometry_dialogs.xenopusBallDialog(currEmbryo.geometry,self).exec_()
			elif "cone" in currEmbryo.geometry.typ:
				ret=pyfrp_gui_geometry_dialogs.coneDialog(currEmbryo.geometry,self).exec_()	
			else:
				ret=pyfrp_gui_geometry_dialogs.geometryDialog(currEmbryo.geometry,self).exec_()	
		return
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Open Geometry Editor for current embryo
		
	def editGeoFile(self):
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			ret=pyfrp_gui_gmsh_editor.gmshFileEditor(currEmbryo.geometry,self).exec_()
		
		return
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Draw Geometry in plot tab
	
	def drawGeometry(self,ann=False):
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			self.createPlotTab('xyz',plotName='Geometry',proj=['3d'])
			currEmbryo.geometry.plotGeometry(ax=self.ax,ann=ann)
			self.adjustCanvas()
		return 
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Analysis handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Edit Analysis Settings
		
	def editAnalysis(self):
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			if currEmbryo.analysis==None:
				currEmbryo.newAnalysis()
			
			ret=pyfrp_gui_analysis_dialogs.analysisDialog(currEmbryo.analysis,self).exec_()
		
		return
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Analyze embryo task/progressbar
	
	def analyzeEmbryo(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		
		
		#Launching Dialog
		self.editAnalysis()
		
		#Genereate wait popup
		self.progressDialog=pyfrp_gui_analysis_dialogs.analysisProgressDialog(None)
		
		#Make backup copy of embryo
		self.originalObj=currEmbryo
		self.backupObj=cpy.deepcopy(currEmbryo)
		
		self.statusBar().showMessage("Analyzing Dataset " + currEmbryo.name)
		
		#Generate Qthread and pass analysis there
		self.task=pyfrp_gui_basics.pyfrpThread()
		self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.analysis.run,signal=self.task.progressSignal)
		
		#Init and start
		self.initTask()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Plots data analysis results of all ROIs of embryo
	
	def plotAllDataTSOfEmbryo(self):
		
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo==None:
			return
		
		self.createPlotTab("intensityTS",plotName=currEmbryo.name+" data",size=[1,1])
		currEmbryo.plotAllData(ax=self.ax)
		
		self.adjustCanvas()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Simulation handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Edit Simulation Settings
		
	def editSimulation(self):
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			if currEmbryo.simulation==None:
				currEmbryo.newSimulation()
			
			ret=pyfrp_gui_simulation_dialogs.simulationSettingsDialog(currEmbryo.simulation,self).exec_()
		
		return
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Simulate embryo task/progressbar
	
	def simulateEmbryo(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
	
		#Launching Dialog
		self.editSimulation()
		
		#Genereate wait popup
		self.progressDialog=pyfrp_gui_simulation_dialogs.simulationProgressDialog(None)
		
		#Make backup copy of embryo
		self.originalObj=currEmbryo
		self.backupObj=cpy.deepcopy(currEmbryo)
		
		self.statusBar().showMessage("Simulating Dataset " + currEmbryo.name)
		
		#Generate Qthread and pass analysis there
		self.task=pyfrp_gui_basics.pyfrpThread()
		self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.simulation.run,signal=self.task.progressSignal)
		
		#Init and start
		self.initTask()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Edit Mesh Settings
		
	def editMesh(self):
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo!=None:
			if currEmbryo.simulation==None:
				currEmbryo.newSimulation()
			ret=pyfrp_gui_mesh_dialogs.meshSettingsDialog(currEmbryo.simulation.mesh,self).exec_()
		return
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Generate Mesh task/progressbar
	
	def generateMesh(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
	
		#Launching Dialog
		self.editMesh()
		
		#Genereate wait popup
		self.progressDialog=pyfrp_gui_mesh_dialogs.genMeshProgressDialog(None)
		
		#Make backup copy of embryo
		self.originalObj=currEmbryo
		self.backupObj=cpy.deepcopy(currEmbryo)
		
		self.statusBar().showMessage("Generating Mesh " + currEmbryo.name)
		
		#Generate Qthread and pass analysis there
		self.task=pyfrp_gui_basics.pyfrpThread()
		self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.simulation.mesh.genMesh)
		
		#Print out mesh stats
		self.printMeshStats()
		
		#Init and start
		self.initTask()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Refine Mesh task/progressbar
	
	def refineMesh(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
	
		#Genereate wait popup
		self.progressDialog=pyfrp_gui_mesh_dialogs.refineMeshProgressDialog(None)
		
		#Make backup copy of embryo
		self.originalObj=currEmbryo
		self.backupObj=cpy.deepcopy(currEmbryo)
		
		self.statusBar().showMessage("Refining Mesh " + currEmbryo.name)
	
		#Generate Qthread and pass analysis there
		self.task=pyfrp_gui_basics.pyfrpThread()
		self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.simulation.mesh.refine)
		
		#Init and start
		self.initTask()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#force Mesh task/progressbar
	
	def forceMeshDensity(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		#Get options
		forceMeshDialog=pyfrp_gui_mesh_dialogs.forceMeshSettingsDialog(currEmbryo.simulation.mesh,self)
		if forceMeshDialog.exec_():
			roiUsed,density,stepPercentage,debug,findIdxs,method,maxCells=forceMeshDialog.getVals()
		
		#Genereate wait popup
		self.progressDialog=pyfrp_gui_mesh_dialogs.forceMeshProgressDialog(None)
		
		#Make backup copy of embryo
		self.originalObj=currEmbryo
		self.backupObj=cpy.deepcopy(currEmbryo)
		
		self.statusBar().showMessage("Computing new mesh with required density " + currEmbryo.name)
		
		#Generate Qthread and pass analysis there
		self.task=pyfrp_gui_basics.pyfrpThread()
		self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.simulation.mesh.forceMinMeshDensityInROI,roiUsed,density,stepPercentage,debug=debug,findIdxs=findIdxs,method=method,maxCells=maxCells)
		
		#Init and start
		self.initTask()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#force Mesh in ROI
	
	def forceROIMeshDensity(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		#Get options
		ret=pyfrp_gui_mesh_dialogs.refineROIMeshSettingsDialog(currEmbryo.simulation.mesh,self).exec_()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Refine Mesh task/progressbar
	
	def printMeshStats(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		currEmbryo.simulation.mesh.printStats()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Plot Mesh
	
	def plotMesh(self):
		
		#Grab embryo 
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		self.createVtkTab()
		
		self.renderer=currEmbryo.simulation.mesh.importVTKFile(sub=True)
		
		self.vtkCanvas.GetRenderWindow().AddRenderer(self.renderer)
		self.vtkCanvas.GetRenderWindow().GetInteractor().Initialize()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Plot Mesh Density
	
	def plotMeshDensity(self):
		
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		self.createPlotTab("meshDensity",plotName=currEmbryo.name+" mesh density",size=[2,2])
		
		currEmbryo.simulation.mesh.plotDensity(axes=self.axes)
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Plots simulation results of all ROIs of embryo
	
	def plotAllSimTSOfEmbryo(self):
		
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo==None:
			return
		
		self.createPlotTab("intensityTS",plotName=currEmbryo.name+" data",size=[1,1])
		currEmbryo.plotAllSim(ax=self.ax)
		
		self.adjustCanvas()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Plots simulation results of all ROIs of embryo
	
	def plotAllSimAndDataTSOfEmbryo(self):
		
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo==None:
			return
		
		self.createPlotTab("intensityTS",plotName=currEmbryo.name+" data",size=[1,1])
		currEmbryo.plotAllData(ax=self.ax)
		currEmbryo.plotAllSim(ax=self.ax)
		
		self.adjustCanvas()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Pinning handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Default Pin
	
	def defaultPinEmbryo(self):
		
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo==None:
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Launch dialog
		ret=pyfrp_gui_pinning_dialogs.defaultPinningDialog(currEmbryo,self).exec_()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Ideal Pin
	
	def idealPinEmbryo(self):
		
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo==None:
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Launch dialog
		ret=pyfrp_gui_pinning_dialogs.idealPinningDialog(currEmbryo,self).exec_()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Fit handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Creates new fit
	
	def newFit(self):
		
		currEmbryo=self.getCurrentEmbryo()
		
		if currEmbryo==None:
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Generate name
		fitNames=pyfrp_misc_module.objAttrToList(currEmbryo.fits,"name")
		newName=pyfrp_misc_module.enumeratedName("newFit",fitNames,sep="_")
		
		#Generate new fit
		newFit=currEmbryo.newFit(newName)
		
		#Launch editor
		ret=pyfrp_gui_fit_dialogs.fitSettingsDialog(newFit,self).exec_()
		
		#Update Object Bar
		newNode=self.updateFitsNodeChildren()
		self.objectBar.setCurrentItem(newNode)
		
		#Perform the new fit
		self.performFit()
		
		return newFit
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Opens edit fit dialog
	
	def editFit(self):
		
		self.getCurrentObj()
		
		if self.currNodeType!='fit':
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Launch editor
		ret=pyfrp_gui_fit_dialogs.fitSettingsDialog(self.currObj,self).exec_()
		
		#Update Object Bar
		self.updateFitsNodeChildren()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Removes fit
	
	def removeFit(self):
		
		self.getCurrentObj()
		
		if self.currNodeType!='fit':
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		currEmbryo=self.getCurrentEmbryo()
		
		currEmbryo.deleteFit(currEmbryo.fits.index(self.currObj))
		
		#Update Object Bar
		self.updateFitsNodeChildren()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Fitting task/progressbar
	
	def performFit(self):
		
		self.getCurrentObj()
		
		if self.currNodeType!='fit':
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		#Check if pinned
		if self.currObj.fitPinned:
			if not self.currObj.checkPinned():
				self.idealPinEmbryo()
		
		#Genereate wait popup
		self.progressDialog=pyfrp_gui_fit_dialogs.fittingProgressDialog(None)
		
		#Make backup copy of embryo
		self.originalObj=self.currObj
		self.backupObj=cpy.deepcopy(self.currObj)
		
		self.statusBar().showMessage("Performing fit " + self.currObj.name)
		
		#Generate Qthread and pass analysis there
		self.task=pyfrp_gui_basics.pyfrpThread()
		self.worker=pyfrp_gui_basics.pyfrpWorker(self.currObj.run)
		
		#Init and start
		self.initTask()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Print fit results
	
	def printFitResults(self):
	
		
		self.getCurrentObj()
		
		if self.currNodeType!='fit':
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		self.currObj.printResults()
		
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Plot Fit
	
	def plotFit(self):
		
		self.getCurrentObj()
		
		if self.currNodeType!='fit':
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if self.currObj.isFitted():
		
			self.createPlotTab("intensityTS",plotName=self.currObj.name+" ",size=[1,1])
			self.currObj.plotFit(ax=self.ax)
			
		self.adjustCanvas()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Statistics handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Opens fit selector
	
	def selectFits(self):
	
		self.getCurrentObj()
	
		if self.currMoleculeNode==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		if len(self.currMolecule.embryos)==0:
			QtGui.QMessageBox.critical(None, "Error","Molecule does not have any embryos to average.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
			
		#Open Fit selector dialog 
		ret=pyfrp_gui_statistics_dialogs.fitSelector(self.currMolecule,True,self).exec_()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Opens crucial parameter selector
	
	def selectCrucialParameters(self):
	
		self.getCurrentObj()
	
		if self.currMoleculeNode==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		if len(self.currMolecule.embryos)==0:
			QtGui.QMessageBox.critical(None, "Error","Molecule does not have any embryos to average.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
			
		#Open Fit selector dialog 
		ret=pyfrp_gui_statistics_dialogs.crucialParameterSelector(self.currMolecule,self).exec_()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Opens crucial parameter selector
	
	def summarizeMolecule(self):
		
		self.getCurrentObj()
		
		if self.currMoleculeNode==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		if len(self.currMolecule.selFits)==0:
			self.selectFits()
		
		self.selectCrucialParameters()
		
		self.currMolecule.sumUpResults()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Task handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Initializes task by connecting task and progressDialog to signals
	
	def initTask(self):
		
		#Set Main GUI set disabled
		self.setDisabled(True)
		
		#Connect signals
		self.worker.taskFinished.connect(self.taskFinished)
		self.task.progressSignal.connect(self.updateProgressDialog)
		self.progressDialog.accepted.connect(self.taskCanceled)
		
		#Move worker and start
		self.worker.moveToThread(self.task)
		self.worker.start.emit()
		
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Updates ProgressBar
	
	def updateProgressDialog(self,n):
		self.progressDialog.progressbar.setValue(n)
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Launched when task is finished. Sets Main GUI available again, closes all dialogs, updates ObjectBar
	
	def taskFinished(self):
		
		#Quit thread
		self.task.quit()
		
		#Close dialog and make GUI available again
		self.progressDialog.close()
		self.statusBar().showMessage("Idle")
		self.setEnabled(True)
		
		#Updating analyzed column in ObjectBar
		currNode=self.getCurrentEmbryoNode()
		self.updateEmbryoNodeProps(currNode)
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Launched when task is canceled. Sets Main GUI available again, closes all dialogs, updates ObjectBar
	
	def taskCanceled(self):
		self.statusBar().showMessage("Idle")
		self.setEnabled(True)
		
		self.task.terminate()
		
		#Map back backup
		self.originalObj=cpy.deepcopy(self.backupObj)
		self.backupObj=None
		
		self.progressDialog.close()
		
		#Updating analyzed column in ObjectBar
		currNode=self.getCurrentEmbryoNode()
		self.updateEmbryoNodeProps(currNode)
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#PlotTab handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Create plot canvas

	def createCanvas(self,parent,size=[1,1],titles=None,sup=None,proj=None,tight=False):
		
		h=10000/self.dpi
		v=10000/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.currTab)
		
		self.fig, self.axes = pyfrp_plot_module.makeSubplot(size,titles=titles,tight=tight,sup=sup,proj=proj,fig=self.fig)
	
		self.ax=self.axes[0]
		
		return self.fig,self.canvas,self.ax

	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Adjusts axes spacing and labels

	def adjustAxes(self,plotType):
		
		if plotType=="bar":
			self.fig.subplots_adjust(bottom=0.3)
		
		elif plotType=="xyz":
			for ax in self.axes:
				ax.set_xlabel("x (um)")
				ax.set_ylabel("y (um)")
				ax.set_zlabel("z (um)")
				
		elif plotType=="intensityTS":	
			self.fig.subplots_adjust(right=0.75)
			for ax in self.axes:
				ax.set_xlabel("Time (s)")
				ax.set_ylabel("Intensity (AU)")
		
		elif plotType=="meshDensity":
			c=["x","y","z",""]
			for i,ax in enumerate(self.axes):
				ax.set_xlabel(c[i])
				ax.set_ylabel("cell Volume")
		
		return 
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Create plot tab
	
	def createPlotTab(self,plotType,plotName="",size=[1,1],titles=None,sup=None,proj=None,tight=False):
		
		#Grab embryo (might need to change this if whole molecule plots)
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return
		
		#Create tab name
		tabname=currEmbryo.name+"/"+plotName+"#1"
		
		#Increment tab counter in name
		for i in range(self.plotTabs.count()):
			if self.plotTabs.tabText(i)==tabname:
				nme,nmbr=tabname.split("#")
				nmbr=str(int(nmbr)+1)
				tabname=nme+"#"+nmbr
		
		#Add new Tab
		self.currTab=QtGui.QWidget()	
		self.plotTabs.addTab(self.currTab,tabname)
		
		self.createCanvas(self.currTab,size=size,titles=titles,sup=sup,proj=proj,tight=tight)
		
		self.currTab.typ=plotType
		
		self.adjustAxes(plotType)

		#Append for bookkeeping
		self.tabAxes.append(self.ax)
		self.tabFigs.append(self.fig)
		
		#Update Canvas
		self.adjustCanvas()
		
		#Check for dummy plot tab
		if self.plotTabs.tabText(0)=="PlotTab":
		
			self.plotTabs.removeTab(self.plotTabs.currentIndex())
			self.plotTabs.setTabsClosable(True)
		
		self.plotTabs.setCurrentWidget(self.currTab)
			
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#What happens when tab is changed
	
	def currentTabChanged(self,value):
			
		self.currTab=self.plotTabs.widget(value)
		if len(self.tabAxes)>0:
			
			self.ax=self.tabAxes[value]
			self.fig=self.tabFigs[value]
		
		if self.currTab!=None:
			self.verticalSplitter.refresh()
			self.horizontalSplitter.refresh()
			self.currTab.setHidden(True)
			self.currTab.setVisible(True)
			self.adjustCanvas()	
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#What happens when tab is closed
	
	def currentTabClosed(self,value):
			
		self.currTab=self.plotTabs.widget(value)
		if len(self.tabAxes)>0:
			self.tabAxes.pop(value)
			self.tabFigs.pop(value)
			
		self.plotTabs.removeTab(value)
		
		if self.plotTabs.count()==0:
			self.createDummpyTab()
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Creates dummy tab
	
	def createDummpyTab(self):
		self.currTab=QtGui.QWidget()
		self.firstTab=self.plotTabs.addTab(self.currTab,"PlotTab")
		self.plotTabs.setTabsClosable(False)
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#vtk handling
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Create tab with vtkcanvas inside
	
	def createVtkTab(self,plotName=""):
		
		#Grab embryo (might need to change this if whole molecule plots)
		currEmbryo=self.getCurrentEmbryo()
		if currEmbryo==None:
			return

		#Create tab name
		tabname=currEmbryo.name+"/"+plotName+"#1"
		
		#Increment tab counter in name
		for i in range(self.plotTabs.count()):
			if self.plotTabs.tabText(i)==tabname:
				nme,nmbr=tabname.split("#")
				nmbr=str(int(nmbr)+1)
				tabname=nme+"#"+nmbr
		
		#Add new Tab
		self.currTab=QtGui.QWidget()
		#self.currTab.setStyleSheet("background-color: rgb(255,0,0); margin:5px; border:1px solid rgb(0, 255, 0); ")
		
		self.adjustCanvas()
		self.createVtkCanvas()
		
		#self.tabAxes.append(self.vtkCanvas)
		#self.tabFigs.append(None)
		
		self.plotTabs.addTab(self.currTab,tabname)
		
		#Check for dummy plot tab
		if self.plotTabs.tabText(0)=="PlotTab":
		
			self.plotTabs.removeTab(self.plotTabs.currentIndex())
			self.plotTabs.setTabsClosable(True)
		
		self.plotTabs.setCurrentWidget(self.currTab)
		
		return self.vtkCanvas
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Create vtkcanvas
	
	def createVtkCanvas(self):
		
		self.vtkCanvas = QVTKRenderWindowInteractor(self.currTab)
		self.vtkCanvas.UpdateSize(int(self.currTab.width()),int(self.currTab.height()))
		
		self.vboxVTKCanvas = QtGui.QVBoxLayout()
		self.vboxVTKCanvas.addWidget(self.vtkCanvas,stretch=1)
		
		self.currTab.setLayout(self.vboxVTKCanvas)
		
		return self.vtkCanvas
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Adjust canvas if slider/splitter changes
	
	def adjustCanvas(self):
		
		if hasattr(self,'currTab'):
			h=int(self.horizontalSplitter.sizes()[1])
			v=int(self.verticalSplitter.sizes()[0])	
			self.currTab.resize(h,v)
		
		if hasattr(self,'fig'):
		
			h=float(self.horizontalSplitter.sizes()[1])/float(self.dpi)
			v=float(self.verticalSplitter.sizes()[0])/float(self.dpi)
			
			if hasattr(self.currTab,'currSlider'):
				
				hSlider=float(self.currTab.currSlider.size().width())/float(self.dpi)
				vSlider=float(self.currTab.currSlider.size().height())/float(self.dpi)
				
				v=v-5*vSlider
			
			self.fig.set_size_inches(h,v,forward=False)
			
			self.canvas.draw()
		
		elif hasattr(self,'vtkCanvas'):
			
			h=int(self.horizontalSplitter.sizes()[1])
			v=int(self.verticalSplitter.sizes()[0])
			
			self.vtkCanvas.UpdateSize(h,v)
			
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Make slider plot
	
	#def create_slider_plot_tab(self,plotType):
		
		##Create new tab
		#if plotType in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int"]:
			#tabname=self.curr_mol.name+"/"+self.curr_bkgd.name+"/"+plotType+"#1"
		#elif plotType in ["ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			#tabname=self.curr_embr.name+"/"+plotType+"#1"	
		#elif plotType=="track_fit":
			#tabname=self.curr_embr.name+"/"+self.curr_fit.name+"/"+"track#1"
			
		#else:	
			#tabname="newtab#1"

		#for i in range(self.plot_tabs.count()):
			#if self.plot_tabs.tabText(i)==tabname:
				#nme,nmbr=tabname.split("#")
				#nmbr=str(int(nmbr)+1)
				#tabname=nme+"#"+nmbr
		
		#self.currTab=QtGui.QWidget()	
		#self.plot_tabs.addTab(self.currTab,tabname)
		
		#h=10000/self.dpi
		#v=10000/self.dpi
		#self.fig = Figure( dpi=self.dpi)
		#self.fig.set_size_inches(h,v,forward=True)
		#self.canvas = FigureCanvas(self.fig)
		#self.canvas.setParent(self.currTab)
		
		##self.mpl_toolbar = NavigationToolbar(self.canvas, self.canvas)
		
		#self.currTab.currSlider=QtGui.QSlider(parent=self.currTab)
		#self.currTab.currSlider.setOrientation(Qt.Horizontal)
		
		#self.currTab.typ=plotType
		
		#self.vbox_slider = QtGui.QVBoxLayout()
		#self.vbox_slider.addWidget(self.canvas,stretch=1)
		
		#if plotType in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int"]:
			#self.currTab.currSlider.setRange(0,shape(self.curr_bkgd.bkgd_vals_slice)[0]-1)
			#self.currTab.currSlider.setSingleStep(1)
			#self.connect(self.currTab.currSlider, QtCore.SIGNAL('sliderReleased()'), self.update_slider_bkgd)
			#self.currTab.lbl_time_slider = QtGui.QLabel("img = 0", self)
			#self.currTab.lbl_time_slider.setAlignment(Qt.AlignHCenter)
			
			#self.vbox_slider.addWidget(self.currTab.currSlider,stretch=0)
			#self.vbox_slider.addWidget(self.currTab.lbl_time_slider)
			
		#if plotType in ["ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			#self.currTab.currSlider.setRange(0,shape(self.curr_embr.vals_slice)[0]-1)
			#self.currTab.currSlider.setSingleStep(1)
			#self.connect(self.currTab.currSlider, QtCore.SIGNAL('sliderReleased()'), self.update_slider_bkgd)
			#self.currTab.lbl_time_slider = QtGui.QLabel("img = 0", self)
			#self.currTab.lbl_time_slider.setAlignment(Qt.AlignHCenter)
			
			#self.vbox_slider.addWidget(self.currTab.currSlider,stretch=0)
			#self.vbox_slider.addWidget(self.currTab.lbl_time_slider)	
		
		#elif plotType=="track_fit":
			#self.currTab.currSlider.setRange(0,shape(self.curr_fit.track_squ_fit)[0]-1)
			#self.currTab.currSlider.setSingleStep(1)
			#self.connect(self.currTab.currSlider, QtCore.SIGNAL('valueChanged(int)'), self.update_slider_track)
			#self.currTab.lbl_track_name_D = QtGui.QLabel("D = ", self)
			#self.currTab.lbl_track_name_prod = QtGui.QLabel("prod = ", self)
			#self.currTab.lbl_track_name_degr = QtGui.QLabel("degr = ", self)
			#self.currTab.lbl_track_name_iter = QtGui.QLabel("iter = ", self)
			#self.currTab.lbl_track_D = QtGui.QLabel("", self)
			#self.currTab.lbl_track_prod = QtGui.QLabel("", self)
			#self.currTab.lbl_track_degr = QtGui.QLabel("", self)
			#self.currTab.lbl_track_iter = QtGui.QLabel("", self)
			#self.currTab.grid_track=QtGui.QGridLayout()
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_name_D,0,0)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_D,0,1)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_name_prod,0,2)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_prod,0,3)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_name_degr,0,4)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_degr,0,5)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_name_iter,0,6)
			#self.currTab.grid_track.addWidget(self.currTab.lbl_track_iter,0,7)
			
			#self.vbox_slider.addWidget(self.currTab.currSlider,stretch=0)
			#self.vbox_slider.addLayout(self.currTab.grid_track)
		
		#self.currTab.setLayout(self.vbox_slider)    
		
		#self.ax = self.fig.add_subplot(111)
		#self.tab_axes.append(self.ax)
		#self.tab_figs.append(self.fig)
		#self.adjust_canvas()
		
		#if self.plot_tabs.tabText(0)=="PlotTab":
		
			#self.plot_tabs.removeTab(self.plot_tabs.currentIndex())
			#self.plot_tabs.setTabsClosable(True)
		
		#self.plot_tabs.setCurrentWidget(self.currTab)
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Update plot through slider for bkgd images
	
	#def update_slider_bkgd(self):
		#value=self.currTab.currSlider.value()
		#self.ax.clear()
		#self.ax.draw_artist(self.currTab.img_axes[value])
		#self.currTab.lbl_time_slider.setText("#img = "+str(value+1))
		#self.currTab.setHidden(True)
		#self.currTab.setVisible(True)
			
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Update plot through slider for tracking images
	
	#def update_slider_track(self,value):
		##value=self.currTab.currSlider.value()
		#self.ax.clear()
		
		#if self.curr_fit.fit_pinned==1:
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_pinned_d,'g*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.out_av_data_pinned_d,'r*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.squ_av_data_pinned_d,'b*')
		
		#else:
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_pinned_d,'g*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.out_av_data_pinned_d,'r*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.squ_av_data_pinned_d,'b*')
		
		#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_slice_fit[value],'g--')
		#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_squ_fit[value],'b--')	
		#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_out_fit[value],'r--')
		
		#self.currTab.lbl_track_D.setText(str(self.curr_fit.track_parms[value][0]))
		#self.currTab.lbl_track_prod.setText(str(self.curr_fit.track_parms[value][1]))
		#self.currTab.lbl_track_degr.setText(str(self.curr_fit.track_parms[value][2]))
		
		#self.currTab.lbl_track_iter.setText(str(value))
		#self.canvas.draw()
	
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#Embryo Menubar Methods
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Add fit
	
	#def add_fit(self):
		
		##Check if top level node is a embryo
		#if self.curr_embr_node==None:
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		#else:
			
			#curr_name="fit_0"
			
			##Check if there is already a fit with the name fit_0
			#in_list=1
			
			#while in_list==1:
				##Go through all embryos
				#for i in range(self.curr_fits.childCount()):
					
					#if self.curr_fits.child(i).data(0,0).toString()==curr_name:
						##If exists, split name in nmbr and name and add +1 to number
						#nme,nmbr=curr_name.split("_")
						#nmbr=str(int(nmbr)+1)
						#curr_name=nme+"_"+nmbr
				
				##Actually in list
				#in_list=0
			
			##Add to embryo bar
			#self.curr_fit_node=QtGui.QTreeWidgetItem(self.curr_fits,[curr_name,'','','0'])
			
			##Create fit object
			#self.curr_embr.add_fit(shape(self.curr_embr.fits)[0],curr_name)
			#self.curr_fit=self.curr_embr.fits[-1]
				
			#ret=pyfrp_subwin.fit_dialog(self.curr_fit,self.curr_mol,self.curr_embr,self).exec_()
			
			##Update bkgd name in embryos column
			#self.curr_fit_node.setText(0,self.curr_fit.name)
			
			#self.embryos_list.setCurrentItem(self.curr_fit_node)
			
			##Perform the new fit
			#self.perform_fit()
			
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Copy fit
	
	#def copy_fit(self):
		
		#if self.curr_node.parent()==self.curr_fits:
			
			#curr_name=self.curr_fit.name
			#curr_fit_mumber=self.curr_fit.fit_number
			#self.curr_fit=cpy.deepcopy(self.curr_fit)
			
			
			#if "copy" in curr_name:
				#new_name=curr_name
			#else:
				#new_name=curr_name+"_copy1"
			
			#in_list=1
			#while in_list==1:
				#for fit in self.curr_embr.fits:
					
					#if fit.name==new_name:
						
						#empty,nmbr=new_name.split("_copy")
						#nmbr=str(int(nmbr)+1)
						#new_name=empty+"_copy"+nmbr
						#break
				#in_list=0
			
			#self.curr_fit.name=new_name
			#self.curr_fit.fit_number=self.curr_fit.fit_number+1
			#self.curr_embr.fit_number=self.curr_embr.fit_number+1
			#self.curr_embr.fits.append(self.curr_fit)
			#self.curr_fit_node=QtGui.QTreeWidgetItem(self.curr_fits,[new_name,'','','0'])		
	
	#def copy_fit_for_other_embryo(self,fit,embryo):
		
		#new_fit=cpy.deepcopy(fit)
		#new_fit.embryo=embryo
		
		#for fit in embryo.fits:
			#if fit.name==new_fit.name:
				#new_fit.name=new_fit.name+"_copied_from"+embryo.name
		
		#embryo.fit_number=embryo.fit_number+1
		#new_fit.fit_number=embryo.fit_number-1
		#embryo.fits.append(new_fit)
		
		#for i in range(self.curr_embryos.childCount()):
		
			#if self.curr_embryos.child(i).data(0,0).toString()==embryo.name:
				#curr_embr_node=self.curr_embryos.child(i)
				#curr_fits=curr_embr_node.child(0)
				#QtGui.QTreeWidgetItem(curr_fits,[new_fit.name,'','','0'])
		
	#def copy_fit_to_all(self):
		
		#if self.curr_fit==None:
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#for emb in self.curr_mol.embryos:
			
			#if emb!=self.curr_embr:
			
				#self.copy_fit_for_other_embryo(self.curr_fit,emb)
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Delete fit
	
	#def delete_fit(self):
		

		##Check what is highlighted
		#if self.embryos_list.currentItem().parent().data(0,0).toString()!="Fits":
			
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
									
		##Remove from list of embryo objects		
		#self.curr_embr.fits.remove(self.curr_fit)
		#self.curr_embr.fit_number=self.curr_embr.fit_number-1
		
		
		##Remove from embryo bar
		#self.parent_node.removeChild(self.curr_node)
		#self.prop_list.clear()
		
		#self.curr_node=self.embryos_list.currentItem()
		#if self.curr_node!=None:
			#self.show_embryo_props()
		#else:
			#self.prop_list.clear()
					
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Edit fit
	
	#def edit_fit(self):
		
		##Check if highlighted node is child of embryo
		#if self.embryos_list.currentItem().parent()!=self.curr_fits:
			
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		#else:
			
			##Open fit dialog
			#ret=pyfrp_subwin.fit_dialog(self.curr_fit,self.curr_mol,self.curr_embr,self).exec_()
			
			##Update property column
			#self.show_embryo_props()
			
			##Update fit name in embryos column
			#self.curr_fit_node=self.embryos_list.currentItem()
			#self.curr_fit_node.setText(0,self.curr_fit.name)
			
		
	#def edit_mult_fit(self):
		
		##Check if highlighted node is child of embryo
		#if self.curr_mol==None:
			
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		#else:
			
			##Backup saving sel_fits
			#sel_fits_backup=list(self.curr_mol.sel_fits)
			
			##Open fit selection dialog
			#ret=pyfrp_subwin.select_fits(self.curr_mol,0,self).exec_()
			
			#if self.curr_mol.sel_fits!=[]:
				##Open mult fit editing dialog
				#ret=pyfrp_subwin.mult_fit_dialog(self.curr_mol,self).exec_()
				
			##Mapping back sel_fits_backup
			#self.curr_mol.sel_fits=list(sel_fits_backup)
			
			##Update property column
			#self.show_embryo_props()
			
	
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Edit Molecule
	
	#def edit_molecule(self):
		
		##Check if highlighted node is child of embryo
		#if self.curr_mol_node==None:
			
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		#else:
			##Open molecule dialog
			#ret=pyfrp_subwin.molecule_dialog(self.curr_mol,self).exec_()
			
			##Update mol name in embryos column
			#self.curr_mol_node.setText(0,self.curr_mol.name)		
			
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Copy embryo
	
	#def copy_embryo(self):
		
		##Check if an molecule is selected
		#if self.curr_embr_node!=None:
			
			#name=str(self.curr_embr.name)
			#curr_name=name+"_copy"
			
			#self.curr_embr=cpy.deepcopy(self.curr_embr)
			#self.curr_embr.name=curr_name
			#self.curr_mol.embryos.append(self.curr_embr)
			
			#emb=self.curr_embr
			
			##Add to embryo bar
			#analyzed=str(0)
			#fitted=str(0)
			#simulated=str(0)
			#if shape(emb.out_av_data_d)[0]>1:
				#analyzed=str(1)
			#if shape(emb.out_av_d)[0]>1:
				#simulated=str(1)	
			#if shape(emb.fits)>0:
				#for fit in emb.fits:
					#if shape(fit.fit_av_d)[0]>1:	
						#fitted=str(1)
				
			##Adding embryo to sidebar
			#self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[emb.name,analyzed,simulated,fitted])
			#self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'','',''])
			
			##Add fits if they exists
			#for fit in emb.fits:
				#if shape(fit.fit_av_d)[0]>0:		
					#QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','','1'])
				#else:	
					#QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','','0'])
								
			#self.curr_embr=emb
			#self.embryos_list.expandItem(self.curr_embr_node)
			#self.embryos_list.expandItem(self.curr_fits)
			
		#else:
			#QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return

	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Copy molecule
	
	#def copy_molecule(self):
		
		##Check if an molecule is selected
		#if self.curr_mol_node!=None:
			
			#name=str(self.curr_mol.name)
			#curr_name=name+"_copy"
			
			#self.curr_mol=cpy.deepcopy(self.curr_mol)
			#self.curr_mol.name=curr_name
			#self.molecules.append(self.curr_mol)
			
			##Adding molecule to sidebar
			#self.curr_node=QtGui.QTreeWidgetItem(self.embryos_list,[curr_name,"","",""])
			#self.curr_embryos=QtGui.QTreeWidgetItem(self.curr_node,["Embryos","","",""])
			#self.curr_mol_node=self.curr_node
			
			#for emb in self.curr_mol.embryos:
		
				##Add to embryo bar
				#analyzed=str(0)
				#fitted=str(0)
				#simulated=str(0)
				#if shape(emb.out_av_data_d)[0]>1:
					#analyzed=str(1)
				#if shape(emb.out_av_d)[0]>1:
					#simulated=str(1)	
				#if shape(emb.fits)>0:
					#for fit in emb.fits:
						#if shape(fit.fit_av_d)[0]>1:	
							#fitted=str(1)
					
				##Adding embryo to sidebar
				#self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[emb.name,analyzed,simulated,fitted])
				#self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'','',''])
			
				##Add fits if they exists
				#for fit in emb.fits:
					#if shape(fit.fit_av_d)[0]>0:		
						#QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','','1'])
					#else:	
						#QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','','0'])
				
				
				#self.curr_embr=emb
				#self.embryos_list.expandItem(self.curr_fits)
									
			#self.embryos_list.expandItem(self.curr_embryos)
		
			
		#else:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Edit properties via doubleclick on property
	
	#def edit_prop(self):
		
		##Getting name of property
		#varname,varvalue=self.prop_list.currentItem().text().split("=")
		
		##Getting property
		#for item in vars(self.curr_embr):
			#if varname==item:
				
				##For general purpose, define what the current property is
				#self.curr_prop=vars(self.curr_embr)[str(item)]
				
				##Check what dialogbox to use
				#if isinstance(self.curr_prop,(int)):
					#varvalue, ok=QtGui.QInputDialog.getInt(self, "Set"+varname, varname+"=")
					#varvalue=int(varvalue)
				#elif isinstance(self.curr_prop,(float)):
					#varvalue, ok=QtGui.QInputDialog.getDouble(self, "Set"+varname, varname+"=")
					#varvalue=float(varvalue)
				#elif isinstance(self.curr_prop,(str)):
					#varvalue, ok=QtGui.QInputDialog.getText(self, "Set"+varname, varname+"=")
					#varvalue=str(varvalue)
					
				##Enter the new value into curr_embr
				#vars(self.curr_embr)[str(item)]=varvalue
				#break
		
		#if varname=="name":
			#self.curr_node.setText(0,varvalue)
			
		##Update property bar
		#self.show_embryo_props()
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Analyze complete molecule
	
	#def analyze_all(self):
		
		#if self.curr_mol_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		##Make backup copy of molecule
		#if self.curr_conf.backup_to_file:
			#self.fn_backup=self.lastopen+"/"+self.curr_mol.name+"_backup.pk"
			#self.curr_mol.save_molecule(self.fn_backup)
		#if self.curr_conf.backup_to_mem:
			#self.temp_mol=cpy.deepcopy(self.curr_mol)
		
		
		
		##Genereate wait popup
		#self.wait_popup=pyfrp_subwin.analyze_all_prog(None,molecule=self.curr_mol)
		#self.wait_popup.accepted.connect(self.analyze_all_canceled)
		#self.statusBar().showMessage("Analyzing molecule " + self.curr_mol.name)
		#self.setDisabled(True)
		
		##Generate Qthread and pass analyze there
		#self.analyze_all_task=pyfrp_subwin.analyze_all_thread(molecule=self.curr_mol)
		#self.analyze_all_task.taskFinished.connect(self.analyze_all_finished)
		#self.analyze_all_task.prog_signal.connect(self.analyze_all_print_prog)
		#self.analyze_all_task.start()
	
	#def analyze_all_print_prog(self,n,k):
	
		#self.wait_popup.progressbar.setValue(k*100+n)
	
	#def analyze_all_finished(self):
		
		#self.wait_popup.close()
		#self.statusBar().showMessage("Idle")
		##Setting analyzed=1
		#for i in range(self.curr_embryos.childCount()):
			#self.curr_embryos.child(i).setText(1,"1")
			#for j in range(self.curr_embryos.child(i).childCount()):
				#if self.curr_embryos.child(i).child(j).data(0,0).toString()!="Fits":
					#self.curr_embryos.child(i).child(j).setText(1,"1")
				
		#self.setEnabled(True)	
		
		#if self.curr_conf.backup_to_mem:
			#self.temp_mol=None
		#if self.curr_conf.backup_to_file:
			#os.remove(self.fn_backup)
			
		#return
	
	#def analyze_all_canceled(self):
		#self.statusBar().showMessage("Idle")
		#self.setEnabled(True)
		
		#self.analyze_all_task.terminate()
		#if self.curr_conf.backup_to_file:
			#self.curr_mol=self.curr_mol.load_molecule(self.fn_backup)
			#os.remove(self.fn_backup)
		#if self.curr_conf.backup_to_mem:
			#self.curr_mol=cpy.deepcopy(self.temp_mol)
			#self.temp_mol=None
			
		#self.wait_popup.close()	
		

	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Edit geometry
	
	#def edit_geometry(self):
		
		##Check if highlighted node is child of embryo
		#if self.curr_embr_node==None:
			
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		#else:
			
			##Open Bkgd dialog for bkgd dataset
			#ret=pyfrp_subwin.geometry_dialog(self.curr_embr,self).exec_()
			
			##Update property column
			#self.embryos_list.setCurrentItem(self.curr_embr_node)
			#self.show_embryo_props()
			
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Edit PDE Parameters
	
	#def edit_pde_parms(self):
		
		##Check if highlighted node is child of embryo
		#if self.curr_embr_node==None:
			
			##If not give error msg
			#QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		#else:
			
			##Open Bkgd dialog for bkgd dataset
			#ret=pyfrp_subwin.sim_dialog(self.curr_embr,self).exec_()
			
			##Update property column
			#self.embryos_list.setCurrentItem(self.curr_embr_node)
			#self.show_embryo_props()		
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Simulate embryo
	
	#def simulate_embryo(self):
		
		#if self.curr_embr_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		##Make backup copy of embryo
		#self.backup_emb=cpy.deepcopy(self.curr_embr)
		
		##Genereate wait popup
		#self.wait_popup=pyfrp_subwin.simulation_prog(None)
		#self.wait_popup.accepted.connect(self.simulation_canceled)
		#self.statusBar().showMessage("Simulating Dataset")
		#self.setDisabled(True)
		
		##Generate Qthread and pass simulate there
		#self.simulation_task=pyfrp_subwin.simulation_thread(embryo=self.curr_embr)
		#self.simulation_task.taskFinished.connect(self.simulation_finished)
		#self.simulation_task.prog_signal.connect(self.simulation_print_prog)
		
		#self.simulation_task.start()
	
	#def simulation_print_prog(self,n):
		
		#self.wait_popup.progressbar.setValue(n)
	
	#def simulation_finished(self):
		
		#self.wait_popup.close()
		#self.statusBar().showMessage("Idle")
		##Setting simulated=1
		#self.curr_embr_node.setText(2,"1")
		
		#self.setEnabled(True)
			
		#return
	
	#def simulation_canceled(self):
		#self.statusBar().showMessage("Idle")
		#self.setEnabled(True)
		
		#self.simulation_task.terminate()
		#self.curr_embr=cpy.deepcopy(self.backup_emb)
		#self.backup_emb=None
		#self.wait_popup.close()
	


	#def plot_sim_timeseries(self):
		
		#if self.curr_embr_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#if shape(self.curr_embr.out_av_d)[0]>0:
			#pass
		#else:
			#QtGui.QMessageBox.critical(None, "Error","No simulation results.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#self.create_plot_tab("data")
		
		#self.ax.plot(self.curr_embr.tvec_sim,self.curr_embr.squ_av_d,'b--',label="squ")
		#self.ax.plot(self.curr_embr.tvec_sim,self.curr_embr.out_av_d,'r--',label="out")
		#self.ax.plot(self.curr_embr.tvec_sim,self.curr_embr.slice_av_d,'g--',label="slice")
		##self.ax.ticklabel_format(style='sci', axis='x', scilimits=(0,1))
		#self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		
		#self.adjust_canvas()
	
	#def plot_sim_data_timeseries(self):
		
		#if self.curr_embr_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#if shape(self.curr_embr.out_av_d)[0]>0:
			#pass
		#else:
			#QtGui.QMessageBox.critical(None, "Error","No simulation results.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#self.create_plot_tab("data")
		
		#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.squ_av_data_d,'b-',label="squ_data")
		#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.out_av_data_d,'r-',label="out_data")
		#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-',label="slice_data")
		#self.ax.plot(self.curr_embr.tvec_sim,self.curr_embr.squ_av_d,'b--',label="squ")
		#self.ax.plot(self.curr_embr.tvec_sim,self.curr_embr.out_av_d,'r--',label="out")
		#self.ax.plot(self.curr_embr.tvec_sim,self.curr_embr.slice_av_d,'g--',label="slice")
		
		##self.ax.ticklabel_format(style='sci', axis='x', scilimits=(0,1))
		#self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		
		#self.adjust_canvas()
	
	
	#def plot_embryo_slice_imgs(self):
		#self.plot_embryo_img_series("slice")
		
	#def plot_embryo_ext_imgs(self):
		#self.plot_embryo_img_series("ext")	
	
	#def plot_embryo_int_imgs(self):
		#self.plot_embryo_img_series("int")
	
	#def plot_embryo_masks_embryo(self):
		#self.plot_embryo_img_series("masks_embryo")
		
	#def plot_embryo_masks_ext(self):
		#self.plot_embryo_img_series("masks_ext")	
	
	#def plot_embryo_masks_int(self):
		#self.plot_embryo_img_series("masks_int")
	
	#def plot_embryo_img_series(self,imgtype):
		
		#if self.curr_embr_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		##Check if embryo is already analyzed
		#if shape(self.curr_embr.slice_av_data_d)[0]==0:
			#reply = QtGui.QMessageBox.question(self, 'Message',"Dataset has not been analyzed yet?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
			#if reply == QtGui.QMessageBox.Yes:
				#self.analyze_dataset()
				#return
			#else:
				#return
		#elif shape(self.curr_embr.masks_embryo)[0]==0 and shape(self.curr_embr.slice_av_data_d)[0]>0:
			#reply = QtGui.QMessageBox.question(self, 'Message',"Image data is currently not loaded, do you want to analyze the data again?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
			#if reply == QtGui.QMessageBox.Yes:
				#self.analyze_dataset()
				#return
			#else:
				#return
			
		#self.create_slider_plot_tab(imgtype)
		#self.img_axes=[]
		
		#self.currTab.imgs=[]
		
		#for i in range(shape(self.curr_embr.vals_slice)[0]):
			#if imgtype=="masks_embryo":
				#img=self.curr_embr.masks_embryo[i]
			#elif imgtype=="masks_ext":
				#img=self.curr_embr.masks_ext[i]
			#elif imgtype=="masks_int":
				#img=self.curr_embr.masks_int[i]
			#elif imgtype=="slice":
				#img=self.curr_embr.vals_slice[i]
			#elif imgtype=="ext":
				#img=self.curr_embr.vals_slice[i]*self.curr_embr.masks_ext[i]
			#elif imgtype=="int":
				#img=self.curr_embr.vals_slice[i]*self.curr_embr.masks_int[i]
				
			#self.img_axes.append(self.ax.imshow(img,cmap='jet'))
			#self.currTab.imgs.append(img)
			
		#self.currTab.img_axes=self.img_axes
		
		#self.adjust_canvas()
		#self.ax.draw_artist(self.currTab.img_axes[0])
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Plot fits
	
	#def plot_fit(self):
		
		#if self.curr_node.parent().data(0,0).toString()!="Fits":
			#QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#if shape(self.curr_fit.squ_av_fitted_d)[0]>0:
			#pass
		#else:
			#QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#self.create_plot_tab("fit")
		#if self.curr_fit.fit_pinned==1:
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.out_av_data_pinned_d,'r*',label='out_data')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.out_av_fitted_d,'r--',label='out_fit')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_pinned_d,'g*',label='slice_data')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.slice_av_fitted_d,'g--',label='slice_fit')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.squ_av_data_pinned_d,'b*',label='squ_data')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.squ_av_fitted_d,'b--',label='squ_fit')
		#else:
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.out_av_data_d,'r*',label='out_data')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.out_av_fitted_d,'r--',label='out_fit')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g*',label='slice_data')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.slice_av_fitted_d,'g--',label='slice_fit')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.squ_av_data_d,'b*',label='squ_data')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.squ_av_fitted_d,'b--',label='squ_fit')
			
		#self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		
		#self.adjust_canvas()
		
	#def plot_track_fit(self):
		
		#if self.embryos_list.currentItem()!=self.curr_fit_node:
			#QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#if shape(self.curr_fit.squ_av_fitted_d)[0]>0:
			#pass
		#else:
			#QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#if self.curr_fit.save_track==1:
			#self.create_slider_plot_tab("track_fit")
					
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_slice_fit[-1],'g--')
			
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.out_av_data_d,'r*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_out_fit[-1],'r--')
			
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.squ_av_data_d,'b*')
			#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_squ_fit[-1],'b--')	
		
			#self.currTab.lbl_track_D.setText(str(self.curr_fit.track_parms[-1][0]))
			#self.currTab.lbl_track_prod.setText(str(self.curr_fit.track_parms[-1][1]))
			#self.currTab.lbl_track_degr.setText(str(self.curr_fit.track_parms[-1][2]))
			
			#self.currTab.currSlider.setSliderPosition(shape(self.curr_fit.track_squ_fit)[0]-1)
			
			#self.canvas.draw()
		
		#else: 
			#QtGui.QMessageBox.critical(None, "Error","Did not save fitting tracks.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	

	
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Perform current fit
	
	#def perform_fit(self):
		
		##Genereate wait popup
		#self.wait_popup=pyfrp_subwin.fitting_prog(None)
		#self.wait_popup.accepted.connect(self.fitting_canceled)
		#self.statusBar().showMessage("Fitting")
		#self.setDisabled(True)
		
		##Generate Qthread and pass fitting there
		#self.fitting_task=pyfrp_subwin.fitting_thread(embryo=self.curr_embr,fit=self.curr_fit,gui=self)
		#self.fitting_task.taskFinished.connect(self.fitting_finished)
		#self.fitting_task.start()
				
	#def fitting_finished(self):
		
		#self.wait_popup.close()
		#self.statusBar().showMessage("Idle")
		##Setting fitted=1
		#self.curr_embr_node.setText(3,"1")
		#self.curr_fit_node.setText(3,"1")
		
		#self.embryos_list.setCurrentItem(self.curr_fit_node)
		#self.curr_node=self.curr_fit_node
		##Plotting
		##if self.curr_fit.save_track==1:
			##self.plot_track_fit()
		##else:
		#self.plot_fit()
		
		#self.setEnabled(True)
	
		#return
	
	#def fitting_canceled(self):
		
		#self.setEnabled(True)
		#self.statusBar().showMessage("Idle")
		#self.fitting_task.terminate()
		
		#self.wait_popup.close()
	
	#def perform_fits_molecule(self):
			
		##Genereate wait popup
		#self.wait_popup=pyfrp_subwin.fitting_prog(None)
		#self.wait_popup.accepted.connect(self.fitting_canceled)
		#self.statusBar().showMessage("Fitting")
		#self.setDisabled(True)
		
		##Generate Qthread and pass fitting there
		#self.fitting_task=pyfrp_subwin.fitting_mol_thread(molecule=self.curr_mol,gui=self)
		#self.fitting_task.taskFinished.connect(self.fitting_all_finished)
		#self.fitting_task.start()
	
	#def fitting_all_finished(self):
		
		#self.wait_popup.close()
		#self.statusBar().showMessage("Idle")
	
		##Setting fitted=1
		#for j in range(self.curr_embryos.childCount()):
			#curr_embr_node=self.curr_embryos.child(j)
			#curr_embr_node.setText(3,"1")
			#curr_fits=curr_embr_node.child(0)
			
			#for i in range(curr_fits.childCount()):
				#curr_fits.child(i).setText(3,"1")
			
		#self.setEnabled(True)
		
		#return	
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Show/Hide Console
	
	#def show_console(self):	
		#self.console.setVisible(True)
		#self.verticalSplitter.refresh()
		#self.horizontalSplitter.refresh()
		#self.adjust_canvas()
		#self.curr_conf.term_hidden=False
		
	#def hide_console(self):
		#self.console.setVisible(False)
		#self.verticalSplitter.refresh()
		#self.horizontalSplitter.refresh()
		#self.adjust_canvas()
		#self.curr_conf.term_hidden=True
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Show/Hide Prop List
	
	#def show_proplist(self):	
		#self.prop_list.setVisible(True)
		#self.verticalSplitter.refresh()
		#self.horizontalSplitter.refresh()
		#self.adjust_canvas()
		#self.curr_conf.prop_hidden=False
		
	#def hide_proplist(self):
		#self.prop_list.setVisible(False)
		#self.verticalSplitter.refresh()
		#self.horizontalSplitter.refresh()
		#self.adjust_canvas()
		#self.curr_conf.prop_hidden=True
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Show/Hide Plotting Tab
	
	#def show_plottab(self):	
		#self.plot_tabs.setVisible(True)
		#self.verticalSplitter.refresh()
		#self.horizontalSplitter.refresh()
		#self.adjust_canvas()
		#self.curr_conf.plot_hidden=False
		
	#def hide_plottab(self):
		#self.plot_tabs.setVisible(False)
		#self.verticalSplitter.refresh()
		#self.horizontalSplitter.refresh()
		#self.curr_conf.plot_hidden=True
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Print current memory usage
	
	#def print_mem_usage(self):
		
		#if self.curr_embr_node!=None:
			
			#if self.curr_fit_node==self.embryos_list.currentItem():
				#pyfrp_misc.print_mem_usage(self.curr_fit)
				#return
			#elif self.curr_noise_node==self.embryos_list.currentItem():
				#pyfrp_misc.print_mem_usage(self.curr_noise)
				#return
			#elif self.curr_pre_node==self.embryos_list.currentItem():
				#pyfrp_misc.print_mem_usage(self.curr_pre)
				#return
			#else:
				#pyfrp_misc.print_mem_usage(self.curr_embr)
				#return
		
		#elif self.curr_bkgd!=None:
		
			#if self.curr_bkgd_pre_node==self.embryos_list.currentItem():
				#pyfrp_misc.print_mem_usage(self.curr_bkgd_pre)
				#return
			#else:
				#pyfrp_misc.print_mem_usage(self.curr_bkgd)
				#return
		
		#else:
			#pyfrp_misc.print_mem_usage(self.curr_mol)
			#return
		
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Export plots/movies	
	
	#def export_plot(self):
		
		##Check if plot is selected
		#ind=self.plot_tabs.currentIndex()
		#if self.plot_tabs.tabText(ind)=="PlotTab":
			#QtGui.QMessageBox.critical(None, "Error","No plot selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		##Building filter
		#file_choices = "IMAGES (*pdf *.png *.tif *eps)"
		
		#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,self.tr(file_choices),"*.png")
		#fn_save=str(fn_save)
		#if fn_save=='':
			#return
		#self.lastopen=os.path.dirname(str(fn_save))
		
		#if  self.currTab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			#ind=self.currTab.currSlider.value()
			#matplotlib.image.imsave(fn_save,self.currTab.imgs[ind])
		#else:
			#self.fig.savefig(fn_save)
	
	#def export_plot_series(self):
		
		#ind=self.plot_tabs.currentIndex()
		#if self.plot_tabs.tabText(ind)=="PlotTab":
			#QtGui.QMessageBox.critical(None, "Error","No plot selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#if self.currTab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			
			#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen)
			#fn_save=str(fn_save)
			#if fn_save=='':
				#return
			#self.lastopen=os.path.dirname(str(fn_save))
			
			#if "." in fn_save:
				#fn_save,ending=fn_save.split(".")
			
			#files=[]
			
			##Make folder to put images in
			#os.mkdir(fn_save)
			
			#fn_file=os.path.basename(fn_save)
			
			#for i in range(shape(self.currTab.img_axes)[0]):
			
				#curr_fn_save=fn_save+"/"+fn_file+str(i)+'.tif'
				#matplotlib.image.imsave(curr_fn_save,self.currTab.imgs[i])
				#files.append(curr_fn_save)
			
		#else:
			#QtGui.QMessageBox.critical(None, "Error","No timeseries selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	
	#def export_movie(self):
		
		#ind=self.plot_tabs.currentIndex()
		#if self.plot_tabs.tabText(ind)=="PlotTab":
			#QtGui.QMessageBox.critical(None, "Error","No plot selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
			
		#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.mpg *.avi","*.mpg")
		#fn_save=str(fn_save)
		#if fn_save=='':
			#return
		#self.lastopen=os.path.dirname(str(fn_save))
		
		#if "." in fn_save:
			#fn_save_temp,ending=fn_save.split(".")
		#else:
			#fn_save_temp=fn_save
			#ending="mpg"
			#fn_save=fn_save+"."+ending
		
		#if ending == "avi" or ending=="mpg":
			#encstr="mencoder 'mf://*_tmp.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=mpeg2video  -o "+ fn_save
		#else:
			#encstr="mencoder 'mf://*_tmp.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=mpeg2video  -o "+ fn_save+'.avi'
		
		#files=[]
		#if self.currTab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int","track_fit"]:	
			#if self.currTab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int"]:	
				#for i in range(shape(self.currTab.img_axes)[0]):
				
					#curr_fn_save=fn_save_temp+str(i)+'_tmp'+'.png'
					#matplotlib.image.imsave(curr_fn_save,self.currTab.imgs[i])
					#files.append(curr_fn_save)
			
			#if self.currTab.typ in ["track_fit"]:
				
				#for i in range(shape(self.curr_fit.track_fit)[0]):
					
					#self.ax.clear()
				
					#if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
							
						#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-')
						#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[i],'g--')
							
					#elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
						
						#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r-')
						#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[i],'r--')
						
					#elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
						
						#self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b-')
						#self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[i],'b--')	
					
					#curr_fn_save=fn_save_temp+str(i)+'_tmp'+'.png'
					#self.fig.savefig(curr_fn_save)
					#files.append(curr_fn_save)
		
			#os.chdir(self.lastopen)
			#os.system(encstr)
				
			#for fn in files:
				#os.remove(fn)

		#else:
			#QtGui.QMessageBox.critical(None, "Error","No timeseries selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Export Embryo object to csv file	
	
	#def export_embryo_csv(self):
		
		#if self.curr_embr_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.cvv","*.csv")
		#fn_save=str(fn_save)
		#if fn_save=='':
			#return
		#self.lastopen=os.path.dirname(str(fn_save))
		
		#pyfrp_misc.write_csv_embryo(fn_save,self.curr_embr)
		
	#def export_molecule_csv(self):
		
		#if self.curr_mol_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.csv","*.csv")
		#fn_save=str(fn_save)
		#if fn_save=='':
			#return
		#self.lastopen=os.path.dirname(str(fn_save))
		
		#pyfrp_misc.write_csv_molecule(fn_save,self.curr_mol)	
	
	#def export_fit_to_csv(self):
		
		#if self.curr_fit_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	
		#if self.curr_fit.k_opt==None: 
			#QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.csv","*.csv")
		#fn_save=str(fn_save)
		#if fn_save=='':
			#return
		
		#if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
							
			#pyfrp_misc.write_csv_timeseries(fn_save,[self.curr_embr.slice_av_data_d,self.curr_fit.fit_av_d],["slice_data","fit"])	
			
			
		#elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
			
			#pyfrp_misc.write_csv_timeseries(fn_save,[self.curr_embr.ext_av_data_d,self.curr_fit.fit_av_d],["ext_data","fit"])	
			
		#elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
			
			#pyfrp_misc.write_csv_timeseries(fn_save,[self.curr_embr.int_av_data_d,self.curr_fit.fit_av_d],["int_data","fit"])	
			
	#def export_errorbar_to_csv(self):
		
		#if self.curr_mol==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	
		#if self.curr_mol.tvec_errors==None: 
			#QtGui.QMessageBox.critical(None, "Error","You haven't made an error bar plot yet. Go to Statistics -> Plotting -> Plot normed fit first.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		#fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.csv","*.csv")
		#fn_save=str(fn_save)
		#if fn_save=='':
			#return
			
		#pyfrp_misc.write_csv_timeseries(fn_save,[self.curr_mol.tvec_avg,self.curr_mol.tvec_errors,self.curr_mol.data_av,self.curr_mol.data_errors,self.curr_mol.fit_av],["tvec_avg","tvec_errors","data_av","data_errors","fit_av"])	
		
		#return
	
	##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	##Statistics
	
	#def sumup_molecule(self):
		
		#if self.curr_mol_node==None:
			#QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
	
		#if shape(self.curr_mol.embryos)[0]==0:
			#QtGui.QMessageBox.critical(None, "Error","Molecule does not have any embryos to average.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			#return
		
		##Open Bkgd dialog for bkgd dataset
		#ret=pyfrp_subwin.select_fits(self.curr_mol,1,self).exec_()
		
		#self.curr_mol.sumup_results()
	
	#def plot_Ds_by_fit(self):
		
		##Get some fits and embryos if not already existent
		#if self.curr_mol.sel_fits==[]:
			#self.sumup_molecule()
		
		#Ds=[]
		#names=[]
		#for fit	in self.curr_mol.sel_fits:
			
			#Ds.append(fit.D_opt_mu)
			#names.append(fit.embryo.name)
		
		#N=shape(self.curr_mol.sel_fits)[0]
		#ind= arange(N)
		#width=0.5
		
		#self.create_plot_tab("bar")
		#self.ax.bar(ind+width, Ds, width, color='m')
		#self.ax.set_xticks(ind+1.5*width)
		#self.ax.set_xticklabels(names,rotation=90)
		#self.adjust_canvas()
		
	#def plot_degrs_by_fit(self):
		
		##Get some fits and embryos if not already existent
		#if self.curr_mol.sel_fits==[]:
			#self.sumup_molecule()
		
		#degrs=[]
		#names=[]
		#for fit	in self.curr_mol.sel_fits:
			
			#degrs.append(fit.degr_opt)
			#names.append(fit.embryo.name)
		
		#N=shape(self.curr_mol.sel_fits)[0]
		#ind= arange(N)
		#width=0.5
		
		#self.create_plot_tab("bar")
		#self.ax.bar(ind+width, degrs, width, color='y')
		#self.ax.set_xticks(ind+1.5*width)
		#self.ax.set_xticklabels(names,rotation=90)
		#self.adjust_canvas()	
		
	#def plot_prods_by_fit(self):
		
		##Get some fits and embryos if not already existent
		#if self.curr_mol.sel_fits==[]:
			#self.sumup_molecule()
		
		#prods=[]
		#names=[]
		#for fit	in self.curr_mol.sel_fits:
			
			#prods.append(fit.prod_opt)
			#names.append(fit.embryo.name)
		
		#N=shape(self.curr_mol.sel_fits)[0]
		#ind= arange(N)
		#width=0.5
		
		#self.create_plot_tab("bar")
		#self.ax.bar(ind+width, prods, width, color='c')
		#self.ax.set_xticks(ind+1.5*width)
		#self.ax.set_xticklabels(names,rotation=90)
		#self.adjust_canvas()		
			
	#def plot_all_by_fit(self):
	
		##Get some fits and embryos if not already existent
		#if self.curr_mol.sel_fits==[]:
			#self.sumup_molecule()
		
		#prods=[]
		#degrs=[]
		#Ds=[]
		#names=[]
		#for fit	in self.curr_mol.sel_fits:
			
			#Ds.append(fit.k_opt)
			#degrs.append(fit.ynaught_opt)
			#prods.append(fit.cnaught_opt)
			#names.append(fit.embryo.name)
		
		#N=shape(self.curr_mol.sel_fits)[0]
		#ind= arange(N)
		#width=0.3
		
		#self.create_plot_tab("bar")
		
		#self.ax.bar(ind-width, Ds, width, color='m',label='D')
		#self.ax2=self.ax.twinx()
		
		#self.ax2.bar(ind, prods, width, color='c',label='prod')
		#self.ax2.bar(ind+width, degrs, width, color='y',label='degr')
		
		
		
		#self.ax.legend(loc=2, borderaxespad=0.)	
		#self.ax2.legend(loc=1, borderaxespad=0.)	
		
		#self.ax.set_xticks(ind+0.5*width)
		#self.ax.set_xticklabels(names,rotation=90)
		
		#self.adjust_canvas()		
		
	#def err_data_fit_plot(self):
		
		##Get some fits and embryos if not already existent
		#if self.curr_mol.sel_fits==[]:
			#self.sumup_molecule()
		
		
		#out_norm=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		#squ_norm=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		#slice_norm=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		
		#out_fit_norm=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		#squ_fit_norm=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		#slice_fit_norm=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		
		#times=zeros((shape(self.curr_mol.sel_fits)[0],self.curr_mol.sel_fits[0].embryo.nframes))
		#last_fit=self.curr_mol.sel_fits[0]
		#i=0
		

		#for fit in self.curr_mol.sel_fits:
			##Check if all datasets have same framerate and nframes
			#if fit.embryo.framerate!=last_fit.embryo.framerate:
				#QtGui.QMessageBox.critical(None, "Error","Embryos do not have the same framerate.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				#return
			#if fit.embryo.nframes!=last_fit.embryo.nframes:
				#QtGui.QMessageBox.critical(None, "Error","Embryos do not have the same number of frames.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				#return
			#if fit.fit_pinned!=last_fit.fit_pinned:
				#QtGui.QMessageBox.critical(None, "Error","Fits do not have same fit_pinned property.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				#return
			
			##Normalize everything between 0 and 1
			#if last_fit.fit_pinned==1:
				#out_norm[i,:]=asarray(fit.embryo.out_av_data_pinned_d)
				#squ_norm[i,:]=asarray(fit.embryo.squ_av_data_pinned_d)
				#slice_norm[i,:]=asarray(fit.embryo.slice_av_data_pinned_d)
				
				
			#elif last_fit.fit_pinned==1:
				#out_norm[i,:]=asarray(fit.embryo.out_av_data_d)
				#squ_norm[i,:]=asarray(fit.embryo.squ_av_data_d)
				#slice_norm[i,:]=asarray(fit.embryo.slice_av_data_d)	
			
			#out_fit_norm[i,:]=asarray(fit.out_av_fitted_d)
			#squ_fit_norm[i,:]=asarray(fit.squ_av_fitted_d)
			#slice_fit_norm[i,:]=asarray(fit.slice_av_fitted_d)	
		
			
			#times[i,:]=fit.embryo.tvec_data
			
			#last_fit=fit
			#i=i+1
				
		
		
		#errors_squ=[]
		#errors_out=[]
		#errors_slice=[]
		
		#avgs_squ=[]
		#avgs_out=[]
		#avgs_slice=[]
		
		#avgs_fit_squ=[]
		#avgs_fit_out=[]
		#avgs_fit_slice=[]
		
		#errors_time=[]
		#avgs_time=[]
		#for i in range(last_fit.embryo.nframes):
			#errors_out.append(std(out_norm[:,i]))
			#errors_squ.append(std(squ_norm[:,i]))
			#errors_slice.append(std(slice_norm[:,i]))
				
			#avgs_out.append(mean(out_norm[:,i]))
			#avgs_squ.append(mean(squ_norm[:,i]))
			#avgs_slice.append(mean(slice_norm[:,i]))
			
			#avgs_fit_out.append(mean(out_fit_norm[:,i]))
			#avgs_fit_squ.append(mean(squ_fit_norm[:,i]))
			#avgs_fit_slice.append(mean(slice_fit_norm[:,i]))
				
			#errors_time.append(std(times[:,i]))
			#avgs_time.append(mean(times[:,i]))
			
		#self.create_plot_tab("err")
		
		#self.curr_mol.tvec_avg=avgs_time
		#self.curr_mol.tvec_errors=errors_time
		#self.curr_mol.squ_fitted_av=avgs_fit_squ
		#self.curr_mol.out_fitted_av=avgs_fit_out
		#self.curr_mol.slice_fitted_av=avgs_fit_slice
		#self.curr_mol.squ_data_av=avgs_squ
		#self.curr_mol.out_data_av=avgs_out
		#self.curr_mol.slice_data_av=avgs_slice
		
		#steps=5
	
		#avgs_time=avgs_time[::steps]
		#avgs_out=avgs_out[::steps]
		#avgs_squ=avgs_squ[::steps]
		#avgs_slice=avgs_slice[::steps]
		
		#errors_out=errors_out[::steps]
		#errors_squ=errors_squ[::steps]
		#errors_slice=errors_slice[::steps]
		
		
		#avgs_fit_out=avgs_fit_out[::steps]
		#avgs_fit_squ=avgs_fit_squ[::steps]
		#avgs_fit_slice=avgs_fit_slice[::steps]
		
	
		#self.ax.errorbar(avgs_time,avgs_out,yerr=errors_out,fmt='ro',label='data_av_ext')
		#self.ax.plot(avgs_time,avgs_fit_out,'r--',label='fit_av_ext')	

		#self.ax.errorbar(avgs_time,avgs_squ,yerr=errors_squ,fmt='bo',label='data_av_int')
		#self.ax.plot(avgs_time,avgs_fit_squ,'b--',label='fit_av_int')

		#self.ax.errorbar(avgs_time,avgs_slice,yerr=errors_slice,fmt='go',label='data_av_slice')
		#self.ax.plot(avgs_time,avgs_fit_slice,'g--',label='fit_av_slice')
			
		#self.ax.autoscale(enable=True, axis='x', tight=True)
		#self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)		
		#self.adjust_canvas()
				
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main
			
def main():
		    
	#Creating application
	#font=QtGui.QFont()
	
	app = QtGui.QApplication(sys.argv)
	font=app.font()
	font.setPointSize(12)
	app.setFont(font)
	
	#Check if stout/sterr should be redirected
	try:
		print sys.argv[1]
		redirect=bool(int(sys.argv[1]))
		print redirect
	except:
		redirect=True
	
	mainWin = pyfrp(redirect=redirect)
	mainWin.show()
	
	sys.exit(app.exec_())
		
if __name__ == '__main__':
	main()


