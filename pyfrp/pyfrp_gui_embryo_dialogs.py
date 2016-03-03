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

#PyQT Dialogs for editing embryo datasets
#(1) embryoDialog

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#QT
from PyQt4 import QtGui, QtCore

#PyFRAP GUI classes
import pyfrp_gui_basics

#PyFRAP modules
from pyfrp_term_module import *
import pyfrp_img_module
import pyfrp_misc_module

#Numpy/Scipy
import numpy as np

#Misc 
import os

#===================================================================================================================================
#Dialog for editing embryo datasets
#===================================================================================================================================

class embryoDialog(pyfrp_gui_basics.basicCanvasDialog):
	
	def __init__(self,embryo,parent):
		
		super(embryoDialog,self).__init__(parent)
		
		#Passing embryo
		self.embryo=embryo
		self.parent=parent
		self.nCharDisplayed=50
		
		#Labels
		self.lblName = QtGui.QLabel("Name:", self)
		
		self.lblDataFT = QtGui.QLabel("Data Filetype:", self)
		self.lblDataEnc = QtGui.QLabel("Data Encoding:", self)
		
		self.lblFnDataFolder = QtGui.QLabel("Data Folder:", self)
	
		
		self.lblDataResPx = QtGui.QLabel("Resolution (px):", self)
		self.lblDataResMu = QtGui.QLabel("Resulution (um):", self)
		
		self.lblSliceDepth = QtGui.QLabel("Imaging Depth (um):", self)
		
		self.lblFrameInterval = QtGui.QLabel("Frame Interval (s):", self)
		self.lblnFrames = QtGui.QLabel("number of frames:", self)
		self.lbltStart = QtGui.QLabel("tStart (s):", self)
		self.lbltEnd = QtGui.QLabel("tEnd (s):", self)
		
		self.lblFnDataFolderValue = QtGui.QLabel("", self)
		
		self.updateDataFolderLbl()
		
		#LineEdits
		self.qleName = QtGui.QLineEdit(self.embryo.name)
		self.qleDataResPx = QtGui.QLineEdit(str(self.embryo.dataResPx))
		self.qleDataResMu = QtGui.QLineEdit(str(self.embryo.dataResMu))
		
		self.qleSliceDepth = QtGui.QLineEdit(str(self.embryo.sliceDepthMu))
		
		self.qleFrameInterval = QtGui.QLineEdit(str(self.embryo.frameInterval))
		self.qlenFrames = QtGui.QLineEdit(str(self.embryo.nFrames))
		self.qletStart = QtGui.QLineEdit(str(self.embryo.tStart))
		self.qletEnd = QtGui.QLineEdit(str(self.embryo.tEnd))
		
		self.doubleValid=QtGui.QDoubleValidator()
		
		self.qleDataResPx.setValidator(self.doubleValid)
		self.qleDataResMu.setValidator(self.doubleValid)
		self.qleSliceDepth.setValidator(self.doubleValid)
		self.qleFrameInterval.setValidator(self.doubleValid)
		self.qlenFrames.setValidator(self.doubleValid)
		self.qletStart.setValidator(self.doubleValid)
		self.qletEnd.setValidator(self.doubleValid)
		
		self.qleName.editingFinished.connect(self.setName)
		self.qleDataResMu.editingFinished.connect(self.setDataResMu)
		self.qleDataResPx.editingFinished.connect(self.setDataResPx)
		self.qleSliceDepth.editingFinished.connect(self.setSliceDepth)
		self.qleFrameInterval.editingFinished.connect(self.setFrameInterval)
		self.qletStart.editingFinished.connect(self.settStart)
		
		self.qletEnd.setReadOnly(True)
		self.qlenFrames.setReadOnly(True)
		
		#Combo
		self.comboDataFt = QtGui.QComboBox(self)
		self.comboDataFt.addItem("tif")
		self.comboDataFt.setCurrentIndex(self.comboDataFt.findText(self.embryo.getDataFT()))
		
		self.comboDataEnc = QtGui.QComboBox(self)
		self.comboDataEnc.addItem("uint8")
		self.comboDataEnc.addItem("uint16")
		self.comboDataEnc.setCurrentIndex(self.comboDataEnc.findText(self.embryo.getDataEnc())) 
		
		self.comboDataFt.activated[str].connect(self.setDataFt)   
		self.comboDataEnc.activated[str].connect(self.setDataEnc)   
	
		#Buttons
		self.btnFnDatafolder=QtGui.QPushButton('Change')
		
		self.btnFnDatafolder.connect(self.btnFnDatafolder, QtCore.SIGNAL('clicked()'), self.setFnDatafolder)
			
		#Layout
		self.dataFolderGrid = QtGui.QGridLayout()
		self.dataFolderGrid.addWidget(self.lblFnDataFolderValue,1,1)
		self.dataFolderGrid.addWidget(self.btnFnDatafolder,1,2)
		
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblName,nRows+1,1)
		self.grid.addWidget(self.lblDataFT,nRows+3,1)
		self.grid.addWidget(self.lblDataEnc,nRows+4,1)
		self.grid.addWidget(self.lblFnDataFolder,nRows+5,1)
		
		self.grid.addWidget(self.lblDataResPx,nRows+8,1)
		self.grid.addWidget(self.lblDataResMu,nRows+9,1)
		
		self.grid.addWidget(self.lblSliceDepth,nRows+10,1)
		
		self.grid.addWidget(self.lblFrameInterval,nRows+12,1)
		self.grid.addWidget(self.lblnFrames,nRows+13,1)
		self.grid.addWidget(self.lbltStart,nRows+14,1)
		self.grid.addWidget(self.lbltEnd,nRows+15,1)
		
		self.grid.addWidget(self.qleName,nRows+1,2)
		self.grid.addWidget(self.comboDataFt,nRows+3,2)
		self.grid.addWidget(self.comboDataEnc,nRows+4,2)
		self.grid.addLayout(self.dataFolderGrid,nRows+5,2)
		
		self.grid.addWidget(self.qleDataResPx,nRows+8,2)
		self.grid.addWidget(self.qleDataResMu,nRows+9,2)
		
		self.grid.addWidget(self.qleSliceDepth,nRows+10,2)
		
		self.grid.addWidget(self.qleFrameInterval,nRows+12,2)
		self.grid.addWidget(self.qlenFrames,nRows+13,2)
		self.grid.addWidget(self.qletStart,nRows+14,2)
		self.grid.addWidget(self.qletEnd,nRows+15,2)
		
		self.showFirstDataImg()
		
		self.show()
		
	def setFnDatafolder(self):
		
		folder = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Data Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.embryo.fnDatafolder=folder
		
		self.parent.lastopen=folder
		
		self.embryo.fnDatafolder=pyfrp_misc_module.slashToFn(self.embryo.fnDatafolder)
	
		self.updateEmbryo()
		self.updateTimeQles()
		
		self.updateDataFolderLbl()
		
		self.showFirstDataImg()
		
	def setFnPreImage(self):
		
		fn = str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.parent.lastopen,))
		if fn=='':
			return
		
		self.parent.lastopen,temp=os.path.split(fn)
		
		self.embryo.fnPreimage=fn
		self.updatePreImageLbl()
		
	def updateDataFolderLbl(self):
		self.lblFnDataFolderValue.setText("..."+self.embryo.fnDatafolder[-self.nCharDisplayed:])
	
	def updatePreImageLbl(self):
		self.lblFnPreimageValue.setText("..."+self.embryo.fnPreimage[-self.nCharDisplayed:])
	
	def updateTimeQles(self):
		self.qlenFrames.setText(str(self.embryo.getNFrames()))
		self.qletStart.setText(str(self.embryo.getTStart()))
		self.qletEnd.setText(str(self.embryo.getTEnd()))
			
	def setDataFt(self,text):
		
		self.embryo.setDataFT(str(text))
		
		self.updateEmbryo()
		self.updateTimeQles()
		
		self.showFirstDataImg()
		
	def setDataEnc(self,text):
	
		self.embryo.setDataEnc("uint16")
		self.showFirstDataImg()
		
	def setName(self):
		self.embryo.setName(str(self.qleName.text()))
		return self.embryo.getName()

	def setFrameInterval(self):
		self.embryo.setFrameInterval(float(str(self.qleFrameInterval.text())))
		self.updateEmbryo()
		self.updateTimeQles()
		return self.embryo.getTvecData()
		
	def settStart(self):
		self.embryo.setTStart(float(str(self.qletStart.text())))
		self.updateEmbryo()
		self.updateTimeQles()
		return self.embryo.getTStart()
		
	def setDataResPx(self):
		self.embryo.setDataResPx(float(str(self.qleDataResPx.text())))
		
		self.ax.set_xlim([0, self.embryo.dataResPx])
		self.ax.set_ylim([0, self.embryo.dataResPx])
		self.canvas.draw()
		return self.embryo.getDataResMu()
	
	def setSliceDepth(self):
		self.embryo.setSliceDepthMu(float(str(self.qleSliceDepth.text())))
	
	def setDataResMu(self):
		self.embryo.setDataResMu(float(str(self.qleDataResMu.text())))
		return self.embryo.getDataResPx()
		
	def updateEmbryo(self):
		self.embryo.updateFileList()
		self.embryo.updateNFrames()
		self.embryo.updateTimeDimensions()
		
	def showFirstDataImg(self):
		
		self.embryo.updateFileList()
		if len(self.embryo.fileList)>0:
		
			fnImg=self.embryo.fnDatafolder+self.embryo.fileList[0]
			img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
			
			self.showImg(img)
	
	