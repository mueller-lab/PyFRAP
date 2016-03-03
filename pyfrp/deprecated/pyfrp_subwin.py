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

#Module containing PyQT classes needed for PyFRAP GUI including:
#1) molecule_dialog: Dialog for editing molecule object
#2) dataset_dialog: Dialog for editing embryo object
#3) analysis_dialog: Dialog for analysis options
#4) geometry_dialog: Dialog for geometry options
#5) sim_dialog: Dialog for simulation options
#6) fit_dialog: Dialog for fitting options
#7) mult_fit_dialog: Dialog for analysis options
#8) about_dialog: Dialog for displaying developer information
#9) analyze_all_prog: Dialog displaying analysis progress of molecule
#10) analyze_all_thread: Thread class for analysis fo whole molecule
#11) analyze_prog: Dialog displaying analysis progress of embryo
#12) analyze_thread: Thread class for analysis
#13) simulation_prog: Dialog displaying fitting progress
#14) simulation_thread: Thread class for simulation
#13) fitting_prog: Dialog displaying fitting progress
#14) fitting_thread: Thread class for fitting
#15) fitting_mol_thread: Thread class if molecule gets fitted
#16) select_fits: Dialog for selecting fits for statistics

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

import sys
from numpy import *
from PyQt4 import QtGui, QtCore
import pyfrp_img_module as pyfrp_img
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_stats_module as pyfrp_stats
import pyfrp_fit_module as pyfrp_fit
from pyfrp_term import *

import time
import os, os.path
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as  NavigationToolbar
from matplotlib.figure import Figure
import skimage.io as skiio
	
#===================================================================================================================================
#Dialog for select/edit fit
#===================================================================================================================================

class analysis_dialog(QtGui.QDialog):
	
	def __init__(self,embryo,parent):
		super(analysis_dialog,self).__init__(parent)
		
		
		self.embryo=embryo
		
		self.setMinimumSize(600,300) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_preimage=QtGui.QPushButton('Change')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		self.btn_preimage.connect(self.btn_preimage, QtCore.SIGNAL('clicked()'), self.set_preimage)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lbl_pre = QtGui.QLabel("Prebleach Image", self)
		self.lbl_pre.setFont(boldfont)
		self.lbl_rim_handling = QtGui.QLabel("Rim Handling", self)
		self.lbl_rim_handling.setFont(boldfont)
		self.lbl_debugging = QtGui.QLabel("Debugging Options", self)
		self.lbl_debugging.setFont(boldfont)
	
		self.lbl_rim = QtGui.QLabel("rim width:", self)
		self.lbl_norm_by_pre = QtGui.QLabel("norm by preimage:", self)
		self.lbl_preimage = QtGui.QLabel("preimage", self)
		self.lbl_img_in_domain = QtGui.QLabel("image in domain:", self)
		self.lbl_add_rim_from_radius = QtGui.QLabel("add rim from outside:", self)
		self.lbl_add_rim = QtGui.QLabel("add rim:", self)
		self.lbl_debug = QtGui.QLabel("debug analysis:", self)
		self.lbl_use_conc_rim = QtGui.QLabel("predefine rim concentration:", self)
		self.lbl_conc_rim = QtGui.QLabel("rim concentration:", self)
		self.lbl_gaussian =  QtGui.QLabel("Gaussian:", self)
		self.lbl_gaussian_sigma =  QtGui.QLabel("Gaussian sigma:", self)
		self.lbl_quad_red =  QtGui.QLabel("Flip quadrant:", self)
		self.lbl_flip_before =  QtGui.QLabel("Flip before process:", self)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_debug = QtGui.QCheckBox('', self)
		self.cb_img_in_domain = QtGui.QCheckBox('', self)
		self.cb_add_rim = QtGui.QCheckBox('', self)
		self.cb_add_rim_from_radius = QtGui.QCheckBox('', self)
		self.cb_conc_rim = QtGui.QCheckBox('', self)
		self.cb_norm_by_pre = QtGui.QCheckBox('', self)
		self.cb_gaussian = QtGui.QCheckBox('', self)
		self.cb_quad_red = QtGui.QCheckBox('', self)
		self.cb_flip_before = QtGui.QCheckBox('', self)
		
		if self.embryo.conc_rim==None:
			self.cb_conc_rim.setCheckState(0)
		else: 
			self.cb_conc_rim.setCheckState(2)
		
		if self.embryo.img_in_domain==0:
			self.cb_img_in_domain.setCheckState(0)
		else: 
			self.cb_img_in_domain.setCheckState(2)
		
		if self.embryo.add_rim_img==0:
			self.cb_add_rim.setCheckState(0)
		else: 
			self.cb_add_rim.setCheckState(2)
		
		if self.embryo.add_rim_from_radius==0:
			self.cb_add_rim_from_radius.setCheckState(0)
		else: 
			self.cb_add_rim_from_radius.setCheckState(2)
		
		if self.embryo.debug_analysis==0:
			self.cb_debug.setCheckState(0)
		else: 
			self.cb_debug.setCheckState(2)
		
		if self.embryo.norm_by_pre==0:
			self.cb_norm_by_pre.setCheckState(0)
		else: 
			self.cb_norm_by_pre.setCheckState(2)
			
		self.cb_gaussian.setCheckState(2*int(self.embryo.gaussian))
		self.cb_quad_red.setCheckState(2*int(self.embryo.quad_red))
		self.cb_flip_before.setCheckState(2*int(self.embryo.flip_before_process))
		
		
		self.connect(self.cb_conc_rim, QtCore.SIGNAL('stateChanged(int)'), self.check_conc_rim)
		self.connect(self.cb_img_in_domain, QtCore.SIGNAL('stateChanged(int)'), self.check_img_in_domain)
		self.connect(self.cb_add_rim, QtCore.SIGNAL('stateChanged(int)'), self.check_add_rim)
		self.connect(self.cb_add_rim_from_radius, QtCore.SIGNAL('stateChanged(int)'), self.check_add_rim_from_radius)
		self.connect(self.cb_debug, QtCore.SIGNAL('stateChanged(int)'), self.check_debug)
		self.connect(self.cb_norm_by_pre, QtCore.SIGNAL('stateChanged(int)'), self.check_norm_by_pre)
		self.connect(self.cb_gaussian, QtCore.SIGNAL('stateChanged(int)'), self.check_gaussian)
		self.connect(self.cb_quad_red, QtCore.SIGNAL('stateChanged(int)'), self.check_quad_red)
		self.connect(self.cb_flip_before, QtCore.SIGNAL('stateChanged(int)'), self.check_flip_before)
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_conc_rim = QtGui.QLineEdit(str(self.embryo.conc_rim))
		self.qle_rim = QtGui.QLineEdit(str(self.embryo.rim))
		
		self.qle_gaussian_sigma = QtGui.QLineEdit(str(self.embryo.gaussian_sigma))
		
		
		#if self.embryo.conc_rim==None:	
			#self.qle_LB_k.setReadOnly(True)
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_conc_rim.setValidator(self.double_valid)
		self.qle_rim.setValidator(self.double_valid)
		self.qle_gaussian_sigma.setValidator(self.double_valid)
		
		self.qle_conc_rim.textChanged[str].connect(self.set_conc_rim)
		self.qle_rim.textChanged[str].connect(self.set_rim)
		self.qle_gaussian_sigma.textChanged[str].connect(self.set_gaussian_sigma)
		
			
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		#Headers
		grid.addWidget(self.lbl_pre,0,1,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_rim_handling,0,3,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_debugging,0,5,1,2,Qt.AlignHCenter)
		
		#First column
		grid.addWidget(self.lbl_norm_by_pre,1,1)
		grid.addWidget(self.lbl_preimage,2,1)
		grid.addWidget(self.lbl_gaussian,3,1)
		grid.addWidget(self.lbl_gaussian_sigma,4,1)
		grid.addWidget(self.lbl_quad_red,5,1)
		grid.addWidget(self.lbl_flip_before,6,1)
		
		#Second column
		grid.addWidget(self.cb_norm_by_pre,1,2)
		grid.addWidget(self.btn_preimage,2,2)
		grid.addWidget(self.cb_gaussian,3,2)
		grid.addWidget(self.qle_gaussian_sigma,4,2)
		grid.addWidget(self.cb_quad_red,5,2)
		grid.addWidget(self.cb_flip_before,6,2)
		
		#Third column
		grid.addWidget(self.lbl_img_in_domain,1,3)
		grid.addWidget(self.lbl_add_rim,2,3)
		grid.addWidget(self.lbl_rim,3,3)
		grid.addWidget(self.lbl_add_rim_from_radius,4,3)
		grid.addWidget(self.lbl_use_conc_rim,5,3)
		grid.addWidget(self.lbl_conc_rim,6,3)
		
		#4th column
		grid.addWidget(self.cb_img_in_domain,1,4)
		grid.addWidget(self.cb_add_rim,2,4)
		grid.addWidget(self.qle_rim,3,4)
		grid.addWidget(self.cb_add_rim_from_radius,4,4)
		grid.addWidget(self.cb_conc_rim,5,4)
		grid.addWidget(self.qle_conc_rim,6,4)
		
		#5th column	
		grid.addWidget(self.lbl_debug,1,5)
		
		#6th column
		grid.addWidget(self.cb_debug,1,6)
		
		#Last Column
		grid.addWidget(self.btn_done,8,6)
		
		grid.setColumnStretch(0,1)
		grid.setColumnStretch(7,1)
	
		#grid.setRowStretch(0,1)
		grid.setRowStretch(9,1)
		
		self.setLayout(grid)
		self.setWindowTitle('Data Analysis Parameters')    
		self.show()
	
	def check_norm_by_pre(self,value):
		
		if value==2:	
			self.embryo.norm_by_pre=1
		else:		
			self.embryo.norm_by_pre=0
			
	def check_img_in_domain(self,value):
		
		if value==2:	
			self.embryo.img_in_domain=1
		else:		
			self.embryo.img_in_domain=0
	
	def check_add_rim(self,value):
		
		if value==2:	
			self.embryo.add_rim_img=1
		else:		
			self.embryo.add_rim_img=0
	
	def check_add_rim_from_radius(self,value):
		
		if value==2:	
			self.embryo.add_rim_from_radius=1
		else:		
			self.embryo.add_rim_from_radius=0
	
	def check_conc_rim(self,value):
		
		if value==2:	
			self.embryo.conc_rim=float(self.qle_conc_rim.text())
		else:		
			self.embryo.conc_rim=None
	
	def check_debug(self,value):
		
		if value==2:	
			self.embryo.debug_analysis=1
		else:		
			self.embryo.debug_analysis=0
	
	def check_gaussian(self,value):
		self.embryo.gaussian=bool(int(value)/2)
	
	def check_quad_red(self,value):
		self.embryo.quad_red=bool(int(value)/2)
		
	def check_flip_before(self,value):
		self.embryo.flip_before_process=bool(int(value)/2)
		
	def set_conc_rim(self,text):
		
		if text=='None':
			self.embryo.conc_rim=None
		else:	
			self.embryo.conc_rim=float(str(text))
		
	def set_rim(self,text):
		
		self.embryo.rim=float(str(text))
		
	def set_gaussian_sigma(self,text):
		if len(str(text))>0:
			self.embryo.gaussian_sigma=float(str(text))
		
	def set_preimage(self):
		
		fn = str(QFileDialog.getOpenFileName(self, "Select Preimage"))
		self.embryo.fn_preimage=fn
		
		#Display image path in lbl
		if len(self.embryo.fn_preimage)>50:
			self.lbl_preimage.setText("..."+self.embryo.fn_preimage[-50:])
		else:
			self.lbl_preimage.setText(self.embryo.fn_preimage)
		
	def done_pressed(self):
			
		self.done(1)
		return self.embryo

#===================================================================================================================================
#Dialog for select/edit geometry
#===================================================================================================================================

class geometry_dialog(QtGui.QDialog):
	def __init__(self,embryo,parent):
		super(geometry_dialog,self).__init__(parent)
			
		self.embryo=embryo
			
		self.dpi = 100
		self.setMinimumSize(1000,500) 
		self.resize(1300,500)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_geometry = QtGui.QLabel("Geometry:", self)
		
		self.lbl_img_in_domain =  QtGui.QLabel("Image in domain:", self)
		self.lbl_fish_inner_radius = QtGui.QLabel("Outer radius (pixels):", self)
		self.lbl_fish_outer_radius = QtGui.QLabel("Inner radius (pixels):", self)
		self.lbl_fish_dist = QtGui.QLabel("Distance between centers (pixels):", self)
		
		self.lbl_cyl_radius = QtGui.QLabel("Radius (pixels):", self)
		self.lbl_cyl_height = QtGui.QLabel("Height (pixels):", self)
		
		self.lbl_frog_radius = QtGui.QLabel("Radius (pixels):", self)
		
		self.lbl_slice_height = QtGui.QLabel("Slice height (pixels):", self)
		
		self.lbl_frog_radius.setVisible(False)
		self.lbl_cyl_radius.setVisible(False)
		self.lbl_cyl_height.setVisible(False)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_geom = QtGui.QComboBox(self)
		self.combo_geom.addItem("Fish")
		self.combo_geom.addItem("Cylinder")
		self.combo_geom.addItem("Frog")
		self.combo_geom.setCurrentIndex(0) 
		
		self.combo_geom.activated[str].connect(self.sel_geom)   
	
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_slice_height = QtGui.QLineEdit(str(self.embryo.slice_height_px[0]))
		
		self.qle_fish_inner_radius = QtGui.QLineEdit(str(self.embryo.fish_inradius_px))
		self.qle_fish_outer_radius = QtGui.QLineEdit(str(self.embryo.fish_outradius_px))
		self.qle_fish_dist = QtGui.QLineEdit(str(self.embryo.fish_dist_px))
		
		self.qle_frog_radius = QtGui.QLineEdit(str(self.embryo.frog_radius_px))
		
		self.qle_cyl_radius = QtGui.QLineEdit(str(self.embryo.cylinder_radius_px))
		self.qle_cyl_height = QtGui.QLineEdit(str(self.embryo.cylinder_height_px))
		
		self.qle_frog_radius.setVisible(False)
		self.qle_cyl_radius.setVisible(False)
		self.qle_cyl_height.setVisible(False)
		
		
		
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_slice_height.setValidator(self.double_valid)
		self.qle_fish_inner_radius.setValidator(self.double_valid)
		self.qle_fish_outer_radius.setValidator(self.double_valid)
		self.qle_fish_dist.setValidator(self.double_valid)
		self.qle_frog_radius.setValidator(self.double_valid)
		self.qle_cyl_height.setValidator(self.double_valid)
		self.qle_cyl_radius.setValidator(self.double_valid)
		
		self.qle_slice_height.editingFinished.connect(self.set_sliceheight)
		self.qle_fish_outer_radius.editingFinished.connect(self.set_fish_outer_radius)
		self.qle_frog_radius.editingFinished.connect(self.set_frog_radius)
		self.qle_cyl_height.editingFinished.connect(self.set_cyl_height)
		self.qle_cyl_radius.editingFinished.connect(self.set_cyl_radius)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_img_in_domain = QtGui.QCheckBox('', self)
		
		if self.embryo.img_in_domain==0:
			self.cb_img_in_domain.setCheckState(0)
			self.qle_fish_dist.setReadOnly(False)
			self.qle_fish_inner_radius.setReadOnly(False)
			self.qle_fish_outer_radius.setReadOnly(False)
			
			self.qle_frog_radius.setReadOnly(False)
			
			self.qle_cyl_radius.setReadOnly(False)
			
		else:	
			self.cb_img_in_domain.setCheckState(2)
			#self.qle_fish_dist.setReadOnly(True)
			#self.qle_fish_inner_radius.setReadOnly(True)
			self.qle_fish_outer_radius.setReadOnly(True)
			
			self.qle_frog_radius.setReadOnly(True)
			
			self.qle_cyl_radius.setReadOnly(True)
		
		self.connect(self.cb_img_in_domain, QtCore.SIGNAL('stateChanged(int)'), self.check_img_in_domain)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Plot frame
		#-------------------------------------------------------------------------------------------------------------------
		
		#self.plot_frame = QtGui.QFrame()
		self.plot_frame = QtGui.QWidget()
		self.plot_frame.setMaximumWidth(1)
		#self.plot_frame = self.create_frame(self.plot_frame)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		#General
		grid.addWidget(self.lbl_geometry,1,1) 
		grid.addWidget(self.lbl_slice_height,2,1)
		grid.addWidget(self.lbl_img_in_domain,3,1)
		
		#Fish
		grid.addWidget(self.lbl_fish_inner_radius,4,1)
		grid.addWidget(self.lbl_fish_outer_radius,5,1)
		grid.addWidget(self.lbl_fish_dist,6,1)
		
		#Frog
		grid.addWidget(self.lbl_frog_radius,4,1)
		
		#Cylinder
		grid.addWidget(self.lbl_cyl_radius,4,1)
		grid.addWidget(self.lbl_cyl_height,5,1)
		
		#General
		grid.addWidget(self.combo_geom,1,2)
		grid.addWidget(self.qle_slice_height,2,2)
		grid.addWidget(self.cb_img_in_domain,3,2)
		
		#Fish
		grid.addWidget(self.qle_fish_inner_radius,4,2)
		grid.addWidget(self.qle_fish_outer_radius,5,2)
		grid.addWidget(self.qle_fish_dist,6,2)
		
		#Frog
		grid.addWidget(self.qle_frog_radius,4,2)
		
		#Cylinder
		grid.addWidget(self.qle_cyl_radius,4,2)
		grid.addWidget(self.qle_cyl_height,5,2)
		
		grid.addWidget(self.btn_done,7,1)
		grid.setColumnMinimumWidth(1,200) 
		
		grid.setRowStretch(9,1)
		grid.setRowStretch(0,1)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Create Canvas
		#-------------------------------------------------------------------------------------------------------------------
		
		self.create_canvas()
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.canvas)
		
		#Add everything to Horizontal Box
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(grid)
		self.hbox.addLayout(self.vbox)
		self.setLayout(self.hbox)    
		
		self.build_fish()
		
		self.setWindowTitle('Edit Geometry')    
		self.show()
	
	def create_frame(self,frame):
		
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		frame.setBackgroundRole(QtGui.QPalette.Light)
		#frame.setAutoFillBackground(Trucentere)        
		frame.setLineWidth (1)
		frame.setFrameShadow (frame.Sunken)
		frame.frameRect().setWidth(100)
		
		return frame
		
	def create_canvas(self):
			
		h=500/self.dpi
		v=500/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.plot_frame)
		
		self.ax = self.fig.add_subplot(111,projection='3d')
		
		self.canvas.draw()
		return 
	
	def sel_geom(self,text):
		
		if text=="Fish":
			self.lbl_fish_inner_radius.setVisible(True)
			self.lbl_fish_outer_radius.setVisible(True)
			self.lbl_fish_dist.setVisible(True)	
			
			self.qle_fish_inner_radius.setVisible(True)
			self.qle_fish_outer_radius.setVisible(True)
			self.qle_fish_dist.setVisible(True)	
			
			self.embryo.geometry="Fish" 
			
			self.plot_fish()
				
		else:
			self.lbl_fish_inner_radius.setVisible(False)
			self.lbl_fish_outer_radius.setVisible(False)
			self.lbl_fish_dist.setVisible(False)
		
			self.qle_fish_inner_radius.setVisible(False)
			self.qle_fish_outer_radius.setVisible(False)
			self.qle_fish_dist.setVisible(False)	
		
		if text=="Frog":
			self.lbl_frog_radius.setVisible(True)
			self.qle_frog_radius.setVisible(True)
			self.embryo.geometry="Frog"
			self.plot_frog()
			
			
		else:
			self.lbl_frog_radius.setVisible(False)
			self.qle_frog_radius.setVisible(False)
			
			
		if text=="Cylinder":
			self.lbl_cyl_radius.setVisible(True)
			self.lbl_cyl_height.setVisible(True)	
			
			self.qle_cyl_radius.setVisible(True)
			self.qle_cyl_height.setVisible(True)	
			
			self.embryo.geometry="Cylinder"
			self.plot_cylinder()
			
		else:
			self.lbl_cyl_radius.setVisible(False)
			self.lbl_cyl_height.setVisible(False)	
			
			self.qle_cyl_radius.setVisible(False)
			self.qle_cyl_height.setVisible(False)	
			
	def check_img_in_domain(self,value):		
		
		if value==2:
		
			self.qle_fish_dist.setReadOnly(False)
			self.qle_fish_inner_radius.setReadOnly(False)
			self.qle_fish_outer_radius.setReadOnly(False)
			
			self.qle_frog_radius.setReadOnly(False)
			self.qle_cyl_radius.setReadOnly(False)
			
			self.embryo.img_in_domain=1
			
		else:	
			
			self.qle_fish_outer_radius.setReadOnly(True)
			
			self.qle_frog_radius.setReadOnly(True)
			
			self.qle_cyl_radius.setReadOnly(True)
			
			self.embryo.img_in_domain=0
			
			self.build_fish()
			self.build_frog()
			self.build_cylinder()
			
	
	def set_sliceheight(self):
		
		text=self.qle_slice_height.text()
		
		self.embryo.slice_height_px[0]=float(str(text))
		
		self.build_fish()
		self.build_frog()
		self.build_cylinder()
		
	def set_fish_outer_radius(self):
		
		text=self.qle_fish_outer_radius.text()
		
		self.embryo.fish_outradius_px=float(str(text))
		self.build_fish()
		
	def set_frog_radius(self):
		
		text=self.qle_frog_radius.text()
		
		self.embryo.frog_radius_px=float(str(text))
		self.build_frog()	
		
	def set_cyl_radius(self):
		
		text=self.qle_cyl_radius.text()
		
		self.embryo.cylinder_radius_px=float(str(text))
		self.build_cylinder()		
			
	def set_cyl_height(self):
		
		text=self.qle_cyl_height.text()
		
		self.embryo.cylinder_height_px=float(str(text))
		self.build_cylinder()		
		
	def build_fish(self):
		
		if self.embryo.img_in_domain==0:
		
			self.embryo.fish_outradius_px=(self.embryo.radius_embr_px**2+self.embryo.slice_height_px[0]**2)/(2*(-self.embryo.slice_height_px[0]))
			self.embryo.fish_inradius_px=self.embryo.fish_outradius_px*1.1
			self.embryo.fish_dist_px=sqrt(self.embryo.fish_inradius_px**2-self.embryo.fish_outradius_px**2) 
				
		else:
			
			self.embryo.fish_inradius_px=self.embryo.fish_outradius_px*1.1
			self.embryo.fish_dist_px=sqrt(self.embryo.fish_inradius_px**2-self.embryo.fish_outradius_px**2) 
		
		#Updating mu coordinates	
		self.embryo.fish_outradius_mu=self.embryo.fish_outradius_px*self.embryo.conv_fact
		self.embryo.fish_in_mu=self.embryo.fish_inradius_px*self.embryo.conv_fact
		self.embryo.fish_dist_mu=self.embryo.fish_dist_px*self.embryo.conv_fact
		
		#Writing into qles
		self.qle_fish_dist.setText(str(self.embryo.fish_dist_px))
		self.qle_fish_inner_radius.setText(str(self.embryo.fish_inradius_px))
		self.qle_fish_outer_radius.setText(str(self.embryo.fish_outradius_px))
		
		if self.combo_geom.currentText()=="Fish":
			self.plot_fish()
		
	def build_frog(self):
		
		if self.embryo.img_in_domain==0:
			
			self.embryo.frog_radius_px=(self.embryo.radius_embr_px**2+self.embryo.slice_height_px[0]**2)/(-2*self.embryo.slice_height_px[0])
		
		#Updating mu coordinates	
		self.embryo.frog_radius_mu=self.embryo.frog_radius_px*self.embryo.conv_fact
		
		#Writing into qles
		self.qle_frog_radius.setText(str(self.embryo.frog_radius_px))
		
		if self.combo_geom.currentText()=="Frog":
			self.plot_frog()
		
	def build_cylinder(self):
			
		if self.embryo.img_in_domain==0:
			
			self.embryo.cylinder_radius_px=self.embryo.radius_embr_px
			self.embryo.cylinder_radius_mu=self.embryo.radius_embr_px*self.embryo.conv_fact
			
		
		#Writing into qles
		self.qle_cyl_radius.setText(str(self.embryo.cylinder_radius_px))
		
		if self.combo_geom.currentText()=="Cylinder":
			self.plot_cylinder()
					
	def plot_fish(self):
		
		res_wire=20
		
		#Define Grid vectors
		v_outer=math.acos(((self.embryo.fish_inradius_px**2-self.embryo.fish_outradius_px**2+self.embryo.fish_dist_px**2)/(2*self.embryo.fish_dist_px)-self.embryo.fish_dist_px)/self.embryo.fish_outradius_px)
		v_inner=math.acos(((self.embryo.fish_inradius_px**2-self.embryo.fish_outradius_px**2+self.embryo.fish_dist_px**2)/(2*self.embryo.fish_dist_px))/self.embryo.fish_inradius_px)
	
		#Parametric vectors
		u=r_[0:2*pi:50j]
		v=r_[v_inner:0:50j]
		#v=r_[0:pi:50j]
	
		x_wire_inner=self.embryo.center_embr_px[0]+self.embryo.fish_inradius_px*outer(cos(u),sin(v))
		y_wire_inner=self.embryo.center_embr_px[1]+self.embryo.fish_inradius_px*outer(sin(u),sin(v))
		z_wire_inner=-(self.embryo.fish_dist_px+self.embryo.fish_outradius_px)+self.embryo.fish_inradius_px*outer(ones(size(u)),cos(v))

		#Parametric vectors
		u=r_[0:2*pi:50j]
		v=r_[0:v_outer:50j]
		#v=r_[0:pi:50j]

		x_wire_outer=self.embryo.center_embr_px[0]+self.embryo.fish_outradius_px*outer(cos(u),sin(v))
		y_wire_outer=self.embryo.center_embr_px[1]+self.embryo.fish_outradius_px*outer(sin(u),sin(v))
		z_wire_outer=-self.embryo.fish_outradius_px+self.embryo.fish_outradius_px*outer(ones(size(u)),cos(v))
	
	
		#Finally plot
		self.ax.cla()
	
		self.ax.plot_wireframe(x_wire_outer,y_wire_outer,z_wire_outer,color='b')
		self.ax.plot_wireframe(x_wire_inner,y_wire_inner,z_wire_inner,color='b')
		
		self.show_img()
		self.adjust_zaxis(z_wire_outer)
		self.canvas.draw()
		
	def plot_frog(self):	
		
		res_wire=20
		
		#Parametric vectors
		u=r_[0:2*pi:50j]
		v=r_[0:pi:50j]
		
		X_wire=self.embryo.center_embr_px[0]+self.embryo.frog_radius_px*outer(cos(u),sin(v))
		Y_wire=self.embryo.center_embr_px[1]+self.embryo.frog_radius_px*outer(sin(u),sin(v))
		Z_wire=-self.embryo.frog_radius_px+self.embryo.frog_radius_px*outer(ones(size(u)),cos(v))
		
		#Finally plot
		self.ax.cla()
	
		self.ax.plot_wireframe(X_wire,Y_wire,Z_wire,color='b')
		self.adjust_zaxis(Z_wire)
		self.show_img()
		
		self.canvas.draw()
		
	def plot_cylinder(self):
		
		res_wire=40
		
		#Define Grid vectors
		x_wire=linspace(-self.embryo.cylinder_radius_px+self.embryo.center_embr_px[0]+0.0001,self.embryo.cylinder_radius_px+self.embryo.center_embr_px[0]-0.0001,res_wire)
		z_wire=linspace(0,-self.embryo.cylinder_height_px,res_wire)
		X_wire, Z_wire=meshgrid(x_wire,z_wire)
		Y_wire=self.embryo.center_embr_px[1]+sqrt(self.embryo.cylinder_radius_px**2-(X_wire-self.embryo.center_embr_px[0])**2)
		Y_wire2=-Y_wire+2*self.embryo.center_embr_px[1]
			
		#Finally plot
		self.ax.cla()
	
		self.ax.plot_wireframe(X_wire,Y_wire,Z_wire)
		self.ax.plot_wireframe(X_wire,Y_wire2,Z_wire)
		
		self.adjust_zaxis(z_wire)
		
		self.show_img()
		
		self.canvas.draw()
	
	def adjust_zaxis(self,z_wire):
		
		if z_wire.min()<self.embryo.slice_height_px[0] and z_wire.max()>self.embryo.slice_height_px[0]:
		
			self.ax.set_zlim([z_wire.min(),z_wire.max()])
		
		elif z_wire.min()>self.embryo.slice_height_px[0]:
			
			self.ax.set_zlim([self.embryo.slice_height_px[0],z_wire.max()])
		
		elif z_wire.max()<self.embryo.slice_height_px[0]:
			
			self.ax.set_zlim([z_wire.min(),self.embryo.slice_height_px[0]])
		
		return
			
		
	def show_img(self):
		
		self.file_list= pyfrp_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
	
		#Check if there is a first image
		if shape(self.file_list)[0]>0:
		
			#Grab first picture
			fn=self.embryo.fn_datafolder+self.file_list[0]
			
			#Check if file exists
			if os.path.isfile(fn)==False:
				print "Warning, no img file in fn_datafolder"
				return 
			
			#Load img
			data_img = skiio.imread(fn).astype(self.embryo.data_enc)
			#data_vals=matplotlib.image.imread(fn)
			data_vals=data_img.real
			data_vals=data_vals.astype('float')
			
			#Make sure img array is 2D
			if len(shape(data_vals))>2:
	
				im_ranges=[]
				for k in range(len(shape(data_vals))):
					im_ranges.append(data_vals[:,:,k].max()-data_vals[:,:,k].min())
				
				ind_max=im_ranges.index(max(im_ranges))
				data_vals=data_vals[:,:,ind_max]
				
			down_sampl=ceil(self.embryo.data_res_px/64.)
			
			data_vals=data_vals[0:self.embryo.data_res_px-1:down_sampl,0:self.embryo.data_res_px-1:down_sampl]
			
			res_new=shape(data_vals)[0]
			
			x_grid=linspace(1,self.embryo.data_res_px,res_new)
			y_grid=linspace(1,self.embryo.data_res_px,res_new)
			
			x_grid, y_grid= meshgrid(x_grid,y_grid)
			
			print shape(x_grid), shape(data_vals)
			
			cflevels=linspace(data_vals.min(),data_vals.max(),10)
			
			print cflevels
			
			#Plot img
			self.ax.contourf(x_grid,y_grid,data_vals,levels=cflevels,offset=self.embryo.slice_height_px[0])
			
	def done_pressed(self):
		
		self.done(1)
		return self.embryo

#===================================================================================================================================
#Dialog for select/edit simulation parameters
#===================================================================================================================================

class sim_dialog(QtGui.QDialog):
	
	def __init__(self,embryo,parent):
		super(sim_dialog,self).__init__(parent)
		
		self.embryo=embryo
				
		self.setMinimumSize(600,300) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_fn_mesh=QtGui.QPushButton('Change')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		self.btn_fn_mesh.connect(self.btn_fn_mesh, QtCore.SIGNAL('clicked()'), self.sel_fn_mesh)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lbl_parms = QtGui.QLabel("PDE Parameters", self)
		self.lbl_parms.setFont(boldfont)
		self.lbl_ics = QtGui.QLabel("Initial conditions", self)
		self.lbl_ics.setFont(boldfont)
		self.lbl_avging = QtGui.QLabel("Averaging", self)
		self.lbl_avging.setFont(boldfont)
		self.lbl_mesh = QtGui.QLabel("Mesh", self)
		self.lbl_mesh.setFont(boldfont)
		self.lbl_debugging = QtGui.QLabel("Debugging", self)
		self.lbl_debugging.setFont(boldfont)
		
		self.lbl_diff = QtGui.QLabel("D=", self)
		self.lbl_prod = QtGui.QLabel("prod=", self)
		self.lbl_degr = QtGui.QLabel("degr=", self)
		self.lbl_steps = QtGui.QLabel("steps=", self)
		self.lbl_tstart = QtGui.QLabel("tstart=", self)
		self.lbl_tend = QtGui.QLabel("tend=", self)
		
		self.lbl_IC_mode = QtGui.QLabel("mode=", self)
		self.lbl_rad_steps = QtGui.QLabel("rad_steps (pixels)=", self)
		
		self.lbl_avg_meth = QtGui.QLabel("averaging method=", self)
		self.lbl_small = QtGui.QLabel("avg_small", self)
		self.lbl_int_meth = QtGui.QLabel("integration method=", self)
		self.lbl_add_rim_sim = QtGui.QLabel("add_rim_sim", self)
		self.lbl_int_steps = QtGui.QLabel("integration steps=", self)
		
		
		self.lbl_volsize = QtGui.QLabel("volSize=", self)
		self.lbl_usemesh = QtGui.QLabel("use mesh=", self)
		self.lbl_name_fn_mesh = QtGui.QLabel("meshfile=", self)
		
		if len(self.embryo.fn_mesh)>50:
			self.lbl_fn_mesh = QtGui.QLabel(str(self.embryo.fn_mesh[-50:]), self)
		else:
			self.lbl_fn_mesh = QtGui.QLabel(str(self.embryo.fn_mesh), self)
			
		self.lbl_debug_all = QtGui.QLabel("debug all", self)
		self.lbl_debug_interp = QtGui.QLabel("debug interpolation", self)
		self.lbl_debug_integr = QtGui.QLabel("debug integration", self)
		self.lbl_debug_out = QtGui.QLabel("debug all (just output)", self)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_small = QtGui.QCheckBox('', self)
		self.cb_add_rim_sim = QtGui.QCheckBox('', self)
		
		self.cb_usemesh = QtGui.QCheckBox('', self)
		
		self.cb_debug_all = QtGui.QCheckBox('', self)
		self.cb_debug_interp = QtGui.QCheckBox('', self)
		self.cb_debug_integr = QtGui.QCheckBox('', self)
		self.cb_debug_out = QtGui.QCheckBox('', self)
		
		self.cb_debug_all.setCheckState(2*self.embryo.debug_simulation)
		self.cb_small.setCheckState(2*self.embryo.debug_simulation)
		self.cb_add_rim_sim.setCheckState(2*self.embryo.debug_simulation)
		self.cb_usemesh.setCheckState(2*self.embryo.debug_simulation)
		
		
		
		self.connect(self.cb_debug_all, QtCore.SIGNAL('stateChanged(int)'), self.check_debug_all)
		self.connect(self.cb_small, QtCore.SIGNAL('stateChanged(int)'), self.check_small)
		self.connect(self.cb_add_rim_sim, QtCore.SIGNAL('stateChanged(int)'), self.check_add_rim_sim)
		self.connect(self.cb_usemesh, QtCore.SIGNAL('stateChanged(int)'), self.check_usemesh)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_diff = QtGui.QLineEdit(str(self.embryo.D))
		self.qle_prod = QtGui.QLineEdit(str(self.embryo.prod))
		self.qle_degr = QtGui.QLineEdit(str(self.embryo.degr))
		
		self.qle_steps = QtGui.QLineEdit(str(self.embryo.steps_sim))
		self.qle_tstart = QtGui.QLineEdit(str(self.embryo.tstart))
		self.qle_tend = QtGui.QLineEdit(str(self.embryo.tend))
	
		self.qle_int_steps = QtGui.QLineEdit(str(self.embryo.int_steps))
		self.qle_rad_steps = QtGui.QLineEdit(str(self.embryo.rad_step_px))
		self.qle_volsize = QtGui.QLineEdit(str(self.embryo.volSize_px))
		
			
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_diff.setValidator(self.double_valid)
		self.qle_prod.setValidator(self.double_valid)
		self.qle_degr.setValidator(self.double_valid)
		self.qle_steps.setValidator(self.double_valid)
		self.qle_rad_steps.setValidator(self.double_valid)
		self.qle_volsize.setValidator(self.double_valid)
		self.qle_int_steps.setValidator(self.double_valid)
		
		self.qle_diff.editingFinished.connect(self.set_diff)
		self.qle_prod.editingFinished.connect(self.set_prod)
		self.qle_degr.editingFinished.connect(self.set_degr)
		self.qle_steps.editingFinished.connect(self.set_steps)
		self.qle_int_steps.editingFinished.connect(self.set_int_steps)
		self.qle_volsize.editingFinished.connect(self.set_volsize)
		self.qle_rad_steps.editingFinished.connect(self.set_rad_steps)
		
		self.qle_tstart.setReadOnly(True)
		self.qle_tend.setReadOnly(True)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_ic = QtGui.QComboBox(self)
		self.combo_ic.addItem("ideal")
		self.combo_ic.addItem("radial")
		self.combo_ic.addItem("interpolation")
		self.combo_ic.activated[str].connect(self.sel_ic)   
		
		self.combo_avg = QtGui.QComboBox(self)
		self.combo_avg.addItem("arithmetic average")
		self.combo_avg.addItem("integration")
		self.combo_avg.activated[str].connect(self.sel_avg)   
		
		self.combo_int = QtGui.QComboBox(self)
		self.combo_int.addItem("slow")
		self.combo_int.addItem("fast")
		self.combo_int.addItem("very fast")
		self.combo_int.activated[str].connect(self.sel_int)   
		
		if self.embryo.apply_data==0:
			self.combo_ic.setCurrentIndex(0)
			self.qle_rad_steps.setVisible(False)
		elif self.embryo.apply_data==1:
			self.combo_ic.setCurrentIndex(1)
			self.qle_rad_steps.setVisible(True)
		elif self.embryo.apply_data==3:
			self.combo_ic.setCurrentIndex(2)
			self.qle_rad_steps.setVisible(False)
		
		self.combo_avg.setCurrentIndex(self.embryo.avg_mode)
		if self.embryo.avg_mode==0:
			self.sel_avg("arithmetic average")
		elif self.embryo.avg_mode==1:
			self.sel_avg("integration")
		
		self.combo_int.setCurrentIndex(self.embryo.integration_method)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		#Heading
		grid.addWidget(self.lbl_parms,0,1,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_ics,0,3,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_avging,0,5,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_mesh,0,7,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_debugging,0,10,1,2,Qt.AlignHCenter)
		
		#First column
		grid.addWidget(self.lbl_diff,1,1)
		grid.addWidget(self.lbl_prod,2,1)
		grid.addWidget(self.lbl_degr,3,1)
		grid.addWidget(self.lbl_steps,4,1)
		grid.addWidget(self.lbl_tstart,5,1)
		grid.addWidget(self.lbl_tend,6,1)
		
		#Second column
		grid.addWidget(self.qle_diff,1,2)
		grid.addWidget(self.qle_prod,2,2)
		grid.addWidget(self.qle_degr,3,2)
		grid.addWidget(self.qle_steps,4,2)
		grid.addWidget(self.qle_tstart,5,2)
		grid.addWidget(self.qle_tend,6,2)
		
		#3rd column
		grid.addWidget(self.lbl_IC_mode,1,3)
		grid.addWidget(self.lbl_rad_steps,2,3)
		
		#4th column
		grid.addWidget(self.combo_ic,1,4)
		grid.addWidget(self.qle_rad_steps,2,4)
		
		#5th column
		grid.addWidget(self.lbl_avg_meth,1,5)
		grid.addWidget(self.lbl_int_meth,2,5)
		grid.addWidget(self.lbl_int_steps,3,5)
		grid.addWidget(self.lbl_add_rim_sim,4,5)
		grid.addWidget(self.lbl_small,5,5)
		
		#6th column
		grid.addWidget(self.combo_avg,1,6)
		grid.addWidget(self.combo_int,2,6)
		grid.addWidget(self.qle_int_steps,3,6)
		grid.addWidget(self.cb_add_rim_sim,4,6)
		grid.addWidget(self.cb_small,5,6)
		
		#7th column
		grid.addWidget(self.lbl_usemesh,1,7)
		grid.addWidget(self.lbl_name_fn_mesh,2,7)
		grid.addWidget(self.lbl_volsize,3,7)
		
		#8th column
		grid.addWidget(self.cb_usemesh,1,8)
		grid.addWidget(self.lbl_fn_mesh,2,8)
		grid.addWidget(self.qle_volsize,3,8)
		
		#9th column
		grid.addWidget(self.btn_fn_mesh,2,9)
		
		#10th column
		grid.addWidget(self.lbl_debug_all,1,10)
		grid.addWidget(self.lbl_debug_interp,2,10)
		grid.addWidget(self.lbl_debug_integr,3,10)
		grid.addWidget(self.lbl_debug_out,4,10)
		
		#11th column
		grid.addWidget(self.cb_debug_all,1,11)
		grid.addWidget(self.cb_debug_interp,2,11)
		grid.addWidget(self.cb_debug_integr,3,11)
		grid.addWidget(self.cb_debug_out,4,11)
		
		grid.addWidget(self.btn_done,8,11)
		
		grid.setColumnStretch(0,1)
		grid.setColumnStretch(12,1)
	
		grid.setRowStretch(0,1)
		grid.setRowStretch(9,1)
		
		self.setLayout(grid)    
		self.setWindowTitle('Edit Simulation Parametes')    
		self.show()
	
	def sel_avg(self,text):
		
		if text=="arithmetic average":	
			
			self.embryo.avg_mode=0
			self.lbl_int_meth.setVisible(False)
			self.combo_int.setVisible(False)
			
			self.lbl_int_steps.setVisible(False)
			self.qle_int_steps.setVisible(False)
		
		elif text=="integration":	
	
			self.embryo.avg_mode=1
			self.lbl_int_meth.setVisible(True)
			self.combo_int.setVisible(True)
			
			self.lbl_int_steps.setVisible(True)
			self.qle_int_steps.setVisible(True)
			
	def sel_int(self,text):

		if text=="slow":
			self.embryo.integration_method=0
			
		if text=="fast":
			self.embryo.integration_method=1
		
		if text=="very fast":
			self.embryo.integration_method=2
	
	def sel_ic(self,text):

		if text=="ideal":
			self.embryo.apply_data=0
			self.lbl_rad_steps.setVisible(False)
			
		if text=="radial":
			self.embryo.apply_data=1
			self.lbl_rad_steps.setVisible(True)
			
		if text=="interpolation":
			self.embryo.apply_data=3
			self.lbl_rad_steps.setVisible(False)
	
	def set_rad_steps(self):
		self.embryo.rad_step_px=int(str(self.qle_rad_steps.text()))
		
	def set_int_steps(self):
		self.embryo.int_steps=int(str(self.qle_int_steps.text()))
		
	def set_volsize(self):
		self.embryo.volSize_px=float(str(self.qle_volsize.text()))
		
	def set_diff(self):
		self.embryo.D=float(str(self.qle_diff.text()))
		
	def set_prod(self):
		self.embryo.prod=float(str(self.qle_prod.text()))
	
	def set_degr(self):
		self.embryo.degr=float(str(self.qle_degr.text()))
	
	def set_steps(self):
		self.embryo.steps_sim=int(str(self.qle_steps.text()))
		self.embryo.tvec_sim=linspace(self.embryo.tstart,self.embryo.tend,self.embryo.steps_sim)
		
	def check_usemesh(self,value):
		
		if value==0:
			self.embryo.usemesh=0
			self.lbl_fn_mesh.setVisible(False)
			self.lbl_name_fn_mesh.setVisible(False)
			self.btn_fn_mesh.setVisible(False)
		if value==2:
			self.embryo.usemesh=2
			self.lbl_fn_mesh.setVisible(True)
			self.lbl_name_fn_mesh.setVisible(True)
			self.btn_fn_mesh.setVisible(True)
			
	def check_add_rim_sim(self,value):
		
		if value==0:
			self.embryo.add_rim_sim=0
			
		if value==2:
			self.embryo.add_rim_sim=1
		
		if self.embryo.add_rim_sim!=self.embryo.add_rim_img:
			print "Warning: add_rim_sim is unequal to add_rim_img, for optimal perfomance try to keep them equal."
			
	def check_small(self,value):
		
		if value==0:
			self.embryo.avg_small=0
			
		if value==2:
			self.embryo.avg_small=1
	
	def check_debug_all(self,value):
		
		if value==0:
			self.embryo.debug_simulation=0
			
		if value==2:
			self.embryo.debug_simulation=1
	
	def sel_fn_mesh(self):
		fn = str(QFileDialog.getOpenFileName(self, "Select Meshfile"))
		self.embryo.fn_mesh=fn
		
		#Display image path in lbl
		if len(self.embryo.fn_mesh)>50:
			self.lbl_fn_mesh.setText("..."+self.embryo.fn_mesh[-50:])
		else:
			self.lbl_fn_mesh.setText(self.embryo.fn_mesh)
					
	def done_pressed(self):
			
		self.done(1)
		return self.embryo
	
#===================================================================================================================================
#Dialog for select/edit fit
#===================================================================================================================================

class fit_dialog(QtGui.QDialog):
	
	def __init__(self,fit,molecule,embryo,parent):
		super(fit_dialog,self).__init__(parent)
		
		self.fit=fit
		self.molecule=molecule
		self.embryo=embryo
		
		self.temp_LB_D=[]
		self.temp_UB_D=[]
		self.temp_LB_degr=[]
		self.temp_UB_degr=[]
		self.temp_LB_prod=[]
		self.temp_UB_prod=[]
		
		
		self.setMinimumSize(600,300) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lbl_opt_parms = QtGui.QLabel("Optimization Parameters", self)
		self.lbl_opt_parms.setFont(boldfont)
		self.lbl_guess = QtGui.QLabel("Initial Guess", self)
		self.lbl_guess.setFont(boldfont)
		self.lbl_bounds = QtGui.QLabel("Variable bounds", self)
		self.lbl_bounds.setFont(boldfont)
		self.lbl_bounded = QtGui.QLabel("Bounded?", self)
		self.lbl_bounded.setFont(boldfont)
		self.lbl_fitting = QtGui.QLabel("Fitting options", self)
		self.lbl_fitting.setFont(boldfont)
		
		self.lbl_name = QtGui.QLabel("Name:", self)
		self.lbl_opt_meth = QtGui.QLabel("opt_meth:", self)
		self.lbl_opt_tol = QtGui.QLabel("opt_tol", self)
		self.lbl_maxfun = QtGui.QLabel("maxfun:", self)
		self.lbl_save_track = QtGui.QLabel("save_track:", self)
		self.lbl_debug_fit = QtGui.QLabel("debug_fit:", self)
		
		self.lbl_x0_D = QtGui.QLabel("x0_D:", self)
		self.lbl_x0_degr = QtGui.QLabel("x0_prod:", self)
		self.lbl_x0_prod = QtGui.QLabel("x0_degr:", self)
		
		self.lbl_LB_D = QtGui.QLabel("LB_D:", self)
		self.lbl_UB_D = QtGui.QLabel("UB_D:", self)
		self.lbl_LB_prod = QtGui.QLabel("LB_prod:", self)
		self.lbl_UB_prod = QtGui.QLabel("UB_prod:", self)
		self.lbl_LB_degr = QtGui.QLabel("LB_degr:", self)
		self.lbl_UB_degr = QtGui.QLabel("UB_degr:", self)
		
		self.lbl_fit_out = QtGui.QLabel("fit_out:", self)
		self.lbl_fit_squ = QtGui.QLabel("fit_squ:", self)
		self.lbl_fit_slice = QtGui.QLabel("fit_slice:", self)
		self.lbl_fit_prod = QtGui.QLabel("fit_prod:", self)
		self.lbl_fit_degr = QtGui.QLabel("fit_degr:", self)
		self.lbl_eq = QtGui.QLabel("equalization:", self)
		self.lbl_pin = QtGui.QLabel("fit_pinned:", self)
		
		self.lbl_cut_off = QtGui.QLabel("Cut off at timepoint:", self)
		self.lbl_cut_off_t = QtGui.QLabel("cut_off_t:", self)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_meth = QtGui.QComboBox(self)
		self.combo_meth.addItem("Constrained Nelder-Mead")
		self.combo_meth.addItem("TNC")
		self.combo_meth.addItem("Nelder-Mead")
		self.combo_meth.addItem("L-BFGS-B")
		self.combo_meth.addItem("SLSQP")
		self.combo_meth.addItem("brute")
		self.combo_meth.addItem("BFGS")
		self.combo_meth.addItem("CG")
		self.combo_meth.activated[str].connect(self.sel_meth)   
	
	
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_name = QtGui.QLineEdit(fit.name)
		self.qle_opt_tol = QtGui.QLineEdit(str(self.fit.opt_tol))
		self.qle_maxfun = QtGui.QLineEdit(str(self.fit.maxfun))
		
		self.qle_x0_D = QtGui.QLineEdit(str(self.fit.x0[0]))
		self.qle_x0_degr = QtGui.QLineEdit(str(self.fit.x0[1]))
		self.qle_x0_prod = QtGui.QLineEdit(str(self.fit.x0[2]))
		self.qle_cut_off_t = QtGui.QLineEdit(str(self.fit.cut_off_t))
		
		if fit.LB_D==None:
			self.qle_LB_D = QtGui.QLineEdit("0.0")
			self.qle_LB_D.setReadOnly(True)
		else:	
			
			self.qle_LB_D = QtGui.QLineEdit(str(self.fit.LB_D))
			
		if fit.UB_D==None:
			self.qle_UB_D = QtGui.QLineEdit("0.0")
			self.qle_UB_D.setReadOnly(True)
		else:
			self.qle_UB_D = QtGui.QLineEdit(str(self.fit.UB_D))
			
		if fit.LB_degr==None:
			self.qle_LB_degr = QtGui.QLineEdit("0.0")
			self.qle_LB_degr.setReadOnly(True)
		else:
			self.qle_LB_degr = QtGui.QLineEdit(str(self.fit.LB_degr))
			
		if fit.UB_degr==None:
			self.qle_UB_degr = QtGui.QLineEdit("0.0")
			self.qle_UB_degr.setReadOnly(True)
		else:
			self.qle_UB_degr = QtGui.QLineEdit(str(self.fit.UB_degr))
			
		if fit.LB_prod==None:
			self.qle_LB_prod = QtGui.QLineEdit("0.0")
			self.qle_LB_prod.setReadOnly(True)
		else:
			self.qle_LB_prod = QtGui.QLineEdit(str(self.fit.LB_prod))
			
		if fit.UB_prod==None:			
			self.qle_UB_prod = QtGui.QLineEdit("0.0")
			self.qle_UB_prod.setReadOnly(True)
		else:
			self.qle_UB_prod = QtGui.QLineEdit(str(self.fit.UB_prod))
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_LB_D.setValidator(self.double_valid)
		self.qle_LB_degr.setValidator(self.double_valid)
		self.qle_LB_prod.setValidator(self.double_valid)
		self.qle_UB_D.setValidator(self.double_valid)
		self.qle_UB_degr.setValidator(self.double_valid)
		self.qle_UB_prod.setValidator(self.double_valid)
		self.qle_cut_off_t.setValidator(self.double_valid)
		
		self.qle_name.textChanged[str].connect(self.set_name)
		self.qle_opt_tol.textChanged[str].connect(self.set_opt_tol)
		self.qle_maxfun.textChanged[str].connect(self.set_maxfun)
		
		self.qle_x0_D.textChanged[str].connect(self.set_x0_D)
		self.qle_x0_degr.textChanged[str].connect(self.set_x0_degr)
		self.qle_x0_prod.textChanged[str].connect(self.set_x0_prod)
		
		self.qle_UB_D.textChanged[str].connect(self.set_UB_D)
		self.qle_UB_degr.textChanged[str].connect(self.set_UB_degr)
		self.qle_UB_prod.textChanged[str].connect(self.set_UB_prod)
		
		self.qle_LB_D.textChanged[str].connect(self.set_LB_D)
		self.qle_LB_degr.textChanged[str].connect(self.set_LB_degr)
		self.qle_LB_prod.textChanged[str].connect(self.set_LB_prod)
		
		self.qle_cut_off_t.textChanged[str].connect(self.set_cut_off_t)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_debug_fit = QtGui.QCheckBox('', self)
		self.cb_save_track = QtGui.QCheckBox('', self)
		
		self.cb_bound_LB_D = QtGui.QCheckBox('', self)
		self.cb_bound_LB_degr = QtGui.QCheckBox('', self)
		self.cb_bound_LB_prod = QtGui.QCheckBox('', self)
		
		self.cb_bound_UB_D = QtGui.QCheckBox('', self)
		self.cb_bound_UB_degr = QtGui.QCheckBox('', self)
		self.cb_bound_UB_prod = QtGui.QCheckBox('', self)
		
		self.cb_fit_out = QtGui.QCheckBox('', self)
		self.cb_fit_squ = QtGui.QCheckBox('', self)
		self.cb_fit_slice = QtGui.QCheckBox('', self)
		
		self.cb_fit_degr = QtGui.QCheckBox('', self)
		self.cb_fit_prod = QtGui.QCheckBox('', self)
		
		self.cb_eq = QtGui.QCheckBox('', self)
		self.cb_pin = QtGui.QCheckBox('', self)
		
		self.cb_cut = QtGui.QCheckBox('', self)
		
		self.cb_debug_fit.setCheckState(2*self.fit.debug_fit)	
		self.cb_save_track.setCheckState(2*self.fit.save_track)	
		
		self.cb_eq.setCheckState(2*self.fit.equ_on)	
		self.cb_pin.setCheckState(2*self.fit.fit_pinned)	
		
		self.cb_cut.setCheckState(2*self.fit.fit_cut_off_t)
		
		if fit.LB_D==None:
			self.cb_bound_LB_D.setCheckState(0)
		else: 
			self.cb_bound_LB_D.setCheckState(2)
		
		if fit.UB_D==None:
			self.cb_bound_UB_D.setCheckState(0)
		else: 
			self.cb_bound_UB_D.setCheckState(2)
		
		if fit.LB_degr==None:
			self.cb_bound_LB_degr.setCheckState(0)
		else: 
			self.cb_bound_LB_degr.setCheckState(2)
		
		if fit.UB_degr==None:
			self.cb_bound_UB_degr.setCheckState(0)
		else: 
			self.cb_bound_UB_degr.setCheckState(2)
		
		if fit.LB_prod==None:
			self.cb_bound_LB_prod.setCheckState(0)
		else: 
			self.cb_bound_LB_prod.setCheckState(2)
		
		if fit.UB_prod==None:
			self.cb_bound_UB_prod.setCheckState(0)
		else: 
			self.cb_bound_UB_prod.setCheckState(2)
		
		if fit.fit_cut_off_t==1:
			self.lbl_cut_off_t.setVisible(True)
			self.qle_cut_off_t.setVisible(True)
		else: 
			self.lbl_cut_off_t.setVisible(False)
			self.qle_cut_off_t.setVisible(False)
			
		
		self.connect(self.cb_bound_LB_D, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_D)
		self.connect(self.cb_bound_LB_degr, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_degr)
		self.connect(self.cb_bound_LB_prod, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_prod)
		self.connect(self.cb_bound_UB_D, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_D)
		self.connect(self.cb_bound_UB_degr, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_degr)
		self.connect(self.cb_bound_UB_prod, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_prod)
		
		self.cb_fit_out.setCheckState(2*self.fit.fit_out)
		self.cb_fit_squ.setCheckState(2*self.fit.fit_squ)
		self.cb_fit_slice.setCheckState(2*self.fit.fit_slice)
		self.cb_fit_degr.setCheckState(2*self.fit.fit_degr)
		self.cb_fit_prod.setCheckState(2*self.fit.fit_prod)
		
		self.connect(self.cb_debug_fit, QtCore.SIGNAL('stateChanged(int)'), self.check_debug_fit)
		self.connect(self.cb_save_track, QtCore.SIGNAL('stateChanged(int)'), self.check_save_track)
		
		
		self.connect(self.cb_fit_out, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_out)
		self.connect(self.cb_fit_squ, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_squ)
		self.connect(self.cb_fit_slice, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_slice)
		self.connect(self.cb_fit_degr, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_degr)
		self.connect(self.cb_fit_prod, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_prod)
		
		self.connect(self.cb_eq, QtCore.SIGNAL('stateChanged(int)'), self.check_equ)
		self.connect(self.cb_pin, QtCore.SIGNAL('stateChanged(int)'), self.check_pin)
		self.connect(self.cb_cut, QtCore.SIGNAL('stateChanged(int)'), self.check_cut)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_opt_parms,0,1,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_guess,0,3,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounds,0,6,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounded,0,9,1,1,Qt.AlignHCenter)
		grid.addWidget(self.lbl_fitting,0,10,1,2,Qt.AlignHCenter)
		
		
		grid.addWidget(self.lbl_name,1,1)
		grid.addWidget(self.lbl_opt_meth,2,1)
		grid.addWidget(self.lbl_opt_tol,3,1)
		grid.addWidget(self.lbl_maxfun,4,1)
		grid.addWidget(self.lbl_debug_fit,5,1)
		grid.addWidget(self.lbl_save_track,6,1)
		
		grid.addWidget(self.qle_name,1,2)
		grid.addWidget(self.combo_meth,2,2)
		grid.addWidget(self.qle_opt_tol,3,2)
		grid.addWidget(self.qle_maxfun,4,2)
		grid.addWidget(self.cb_debug_fit,5,2)
		grid.addWidget(self.cb_save_track,6,2)
		
		grid.addWidget(self.lbl_x0_D,1,3)
		grid.addWidget(self.lbl_x0_degr,2,3)
		grid.addWidget(self.lbl_x0_prod,3,3)
		
		grid.addWidget(self.qle_x0_D,1,4)
		grid.addWidget(self.qle_x0_degr,2,4)
		grid.addWidget(self.qle_x0_prod,3,4)
		
		
		
		grid.addWidget(self.lbl_LB_D,1,6)
		grid.addWidget(self.lbl_LB_degr,3,6)
		grid.addWidget(self.lbl_LB_prod,5,6)
		
		grid.addWidget(self.lbl_UB_D,2,6)
		grid.addWidget(self.lbl_UB_degr,4,6)
		grid.addWidget(self.lbl_UB_prod,6,6)
		
		grid.addWidget(self.qle_LB_D,1,7)
		grid.addWidget(self.qle_LB_degr,3,7)
		grid.addWidget(self.qle_LB_prod,5,7)
		
		grid.addWidget(self.qle_UB_D,2,7)
		grid.addWidget(self.qle_UB_degr,4,7)
		grid.addWidget(self.qle_UB_prod,6,7)
	
		
		grid.addWidget(self.cb_bound_LB_D,1,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_degr,3,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_prod,5,9,Qt.AlignHCenter)
		
		grid.addWidget(self.cb_bound_UB_D,2,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_degr,4,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_prod,6,9,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_fit_out,1,10)
		grid.addWidget(self.lbl_fit_squ,2,10)
		grid.addWidget(self.lbl_fit_slice,3,10)
		grid.addWidget(self.lbl_fit_degr,4,10)
		grid.addWidget(self.lbl_fit_prod,5,10)
		grid.addWidget(self.lbl_eq,6,10)
		grid.addWidget(self.lbl_pin,7,10)
		grid.addWidget(self.lbl_cut_off,8,10)
		grid.addWidget(self.lbl_cut_off_t,9,10)
		
		grid.addWidget(self.cb_fit_out,1,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_squ,2,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_slice,3,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_degr,4,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_prod,5,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_eq,6,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_pin,7,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_cut,8,11,Qt.AlignHCenter)
		grid.addWidget(self.qle_cut_off_t,9,11,Qt.AlignHCenter)
		
		
		grid.addWidget(self.btn_done,10,11)
		
		grid.setColumnStretch(0,1)
		grid.setColumnStretch(12,1)
	
		#grid.setRowStretch(0,1)
		grid.setRowStretch(9,1)
		
		self.setLayout(grid)    
		self.setWindowTitle('Edit Fit')    
		self.show()
	
	def set_name(self,text):
		if text!="":
			self.fit.name=str(text)
	
	def set_opt_tol(self,text):
		if text!="":
			self.fit.opt_tol=float(str(text))
		
	def set_maxfun(self,text):
		if text!="":
			self.fit.maxfun=int(str(text))
		
	def set_x0_D(self,text):
		if text!="":
			self.fit.x0[0]=float(str(text))
		
	def set_x0_degr(self,text):
		if text!="":
			self.fit.x0[1]=float(str(text))
			
	def set_x0_prod(self,text):
		if text!="":
			self.fit.x0[2]=float(str(text))
			
	def set_UB_D(self,text):
		if text!="":
			self.fit.UB_D=float(str(text))
		
	def set_UB_degr(self,text):
		if text!="":
			self.fit.UB_degr=float(str(text))
	
	def set_UB_prod(self,text):
		if text!="":
			self.fit.UB_prod=float(str(text))
			
	def set_LB_D(self,text):
		if text!="":
			self.fit.LB_D=float(str(text))
		
	def set_LB_degr(self,text):
		if text!="":
			self.fit.LB_degr=float(str(text))
	
	def set_LB_prod(self,text):
		if text!="":
			self.fit.LB_prod=float(str(text))
	
	def set_cut_off_t(self,text):
		if text!="":
			self.fit.cut_off_t=float(str(text))
	
	def sel_meth(self,text):
		self.fit.opt_meth=str(text)
		self.update_bounds_after_meth(str(text))
		
	def update_bounds_after_meth(self,text):
	
		if text not in ["TNC","L-BFGS-B","SLSQP","brute"]:
			
			self.temp_LB_D=self.fit.LB_D
			self.temp_UB_D=self.fit.UB_D
			self.temp_LB_degr=self.fit.LB_degr
			self.temp_UB_degr=self.fit.UB_degr
			self.temp_LB_prod=self.fit.LB_prod
			self.temp_UB_prod=self.fit.UB_prod
			
			self.check_bound_LB_D(0)
			self.check_bound_UB_D(0)
			self.check_bound_LB_degr(0)
			self.check_bound_UB_degr(0)
			self.check_bound_LB_prod(0)
			self.check_bound_UB_prod(0)
			
			self.cb_bound_LB_D.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_D.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Unchecked)
			
			
		else:
			if self.temp_LB_D==[] or text=="brute":
				
				self.check_bound_LB_D(2)
				self.check_bound_UB_D(2)
				self.check_bound_LB_degr(2)
				self.check_bound_UB_degr(2)
				self.check_bound_LB_prod(2)
				self.check_bound_UB_prod(2)
				
				self.cb_bound_LB_D.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_D.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Checked)
				
				if text=="brute":
					self.bounds_for_brute()
				
			else: 
				if self.temp_LB_D==None:
					self.check_bound_LB_D(0)
					self.cb_bound_LB_D.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_D(2)
					self.cb_bound_LB_D.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_D==None:
					self.check_bound_UB_D(0)
					self.cb_bound_UB_D.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_D(2)
					self.cb_bound_UB_D.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_LB_degr==None:
					self.check_bound_LB_degr(0)
					self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_degr(2)
					self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_degr==None:
					self.check_bound_UB_degr(0)
					self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_degr(2)
					self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_LB_prod==None:
					self.check_bound_LB_prod(0)
					self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_prod(2)
					self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_prod==None:
					self.check_bound_UB_prod(0)
					self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_degr(2)
					self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Checked)
		
		return
		
	def bounds_for_brute(self):

		if self.fit.LB_D>=self.fit.UB_D:
			print "LBs need to be larger than UBs for D, going to fix this."
			if self.fit.UB_D>0:
				self.fit.UB_D==10*self.fit.LB_D
			else:
				self.fit.UB_D=2*self.fit.embryo.D
			
			self.qle_UB_D.setText(str(self.fit.UB_D))
			
		if self.fit.LB_degr>=self.fit.UB_degr:
			print "LBs need to be larger than UBs for degr, going to fix this."
			self.fit.UB_degr=1.5*self.fit.LB_degr
			self.qle_UB_degr.setText(str(self.fit.UB_degr))
			
		if self.fit.LB_prod>=self.fit.UB_prod:
			print "LBs need to be larger than UBs for prod, going to fix this."
			
			self.qle_UB_prod.setText(str(self.fit.UB_prod))
			self.combo_UB_prod.setCurrentIndex(1)
		return
	
	def check_equ(self, value):
		self.fit.equ_on=int(value/2)
	
	def check_pin(self, value):
		self.fit.fit_pinned=int(value/2)
	
	def check_cut(self, value):
		self.fit.fit_cut_off_t=int(value/2)
		self.qle_cut_off_t.setVisible(bool(int(value/2)))
		self.lbl_cut_off_t.setVisible(bool(int(value/2)))
		
	def check_fit_out(self, value):
		self.fit.fit_out=int(value/2)
		
	def check_fit_slice(self, value):
		self.fit.fit_slice=int(value/2)
			
	def check_fit_squ(self, value):
		self.fit.fit_squ=int(value/2)
				
	def check_fit_degr(self, value):
		self.fit.fit_degr=int(value/2)
	
	def check_fit_prod(self, value):
		self.fit.fit_prod=int(value/2)
	
	def check_debug_fit(self, value):
		self.fit.debug_fit=int(value/2)
		
	def check_save_track(self, value):
		self.fit.save_track=int(value/2)	
	
	def check_bound_LB_D(self,value):
		if value==2:
			self.qle_LB_D.setReadOnly(False)
			self.fit.LB_D=float(self.qle_LB_D.text())	
		elif value==0:	
			self.qle_LB_D.setReadOnly(True)
			self.fit.LB_D=None
			
	def check_bound_UB_D(self,value):
		if value==2:
			self.qle_UB_D.setReadOnly(False)
			self.fit.UB_D=float(self.qle_UB_D.text())	
		elif value==0:	
			self.qle_UB_D.setReadOnly(True)
			self.fit.UB_D=None		
	
	def check_bound_LB_degr(self,value):
		if value==2:
			self.qle_LB_degr.setReadOnly(False)
			self.fit.LB_degr=float(self.qle_LB_degr.text())	
		elif value==0:	
			self.qle_LB_degr.setReadOnly(True)
			self.fit.LB_degr=None
			
	def check_bound_UB_degr(self,value):
		if value==2:
			self.qle_UB_degr.setReadOnly(False)
			self.fit.UB_degr=float(self.qle_UB_degr.text())	
		elif value==0:	
			self.qle_UB_degr.setReadOnly(True)
			self.fit.UB_degr=None	
			
	def check_bound_LB_prod(self,value):
		if value==2:
			self.qle_LB_prod.setReadOnly(False)
			self.fit.LB_prod=float(self.qle_LB_prod.text())	
		elif value==0:	
			self.qle_LB_prod.setReadOnly(True)
			self.fit.LB_prod=None
			
	def check_bound_UB_prod(self,value):
		if value==2:
			self.qle_UB_prod.setReadOnly(False)
			self.fit.UB_prod=float(self.qle_UB_prod.text())	
		elif value==0:	
			self.qle_UB_prod.setReadOnly(True)
			self.fit.UB_prod=None			
			
	def done_pressed(self):
		if self.fit.fit_out==0 and self.fit.fit_squ==0 and self.fit.fit_slice==0:
			print "Warning, you are fitting none of the regions, you need to check at least one!"
		else:	
			self.done(1)
		return 


#===================================================================================================================================
#Dialog to edit multiple fits
#===================================================================================================================================

class mult_fit_dialog(QtGui.QDialog):
	
	
	def __init__(self,molecule,parent):
		super(mult_fit_dialog,self).__init__(parent)
		
	
		self.molecule=molecule

		
		self.temp_LB_D=[]
		self.temp_UB_D=[]
		self.temp_LB_degr=[]
		self.temp_UB_degr=[]
		self.temp_LB_prod=[]
		self.temp_UB_prod=[]
		
		self.setMinimumSize(600,300) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lbl_opt_parms = QtGui.QLabel("Optimization Parameters", self)
		self.lbl_opt_parms.setFont(boldfont)
		self.lbl_guess = QtGui.QLabel("Initial Guess", self)
		self.lbl_guess.setFont(boldfont)
		self.lbl_bounds = QtGui.QLabel("Variable bounds", self)
		self.lbl_bounds.setFont(boldfont)
		self.lbl_bounded = QtGui.QLabel("Bounded?", self)
		self.lbl_bounded.setFont(boldfont)
		self.lbl_fitting = QtGui.QLabel("Fitting options", self)
		self.lbl_fitting.setFont(boldfont)
		
		
		self.lbl_opt_meth = QtGui.QLabel("opt_meth:", self)
		self.lbl_opt_tol = QtGui.QLabel("opt_tol", self)
		self.lbl_maxfun = QtGui.QLabel("maxfun:", self)
		self.lbl_save_track = QtGui.QLabel("save_track:", self)
		self.lbl_debug_fit = QtGui.QLabel("debug_fit:", self)
		
		self.lbl_x0_D = QtGui.QLabel("x0_D:", self)
		self.lbl_x0_degr = QtGui.QLabel("x0_prod:", self)
		self.lbl_x0_prod = QtGui.QLabel("x0_degr:", self)
		
		self.lbl_LB_D = QtGui.QLabel("LB_D:", self)
		self.lbl_UB_D = QtGui.QLabel("UB_D:", self)
		self.lbl_LB_prod = QtGui.QLabel("LB_prod:", self)
		self.lbl_UB_prod = QtGui.QLabel("UB_prod:", self)
		self.lbl_LB_degr = QtGui.QLabel("LB_degr:", self)
		self.lbl_UB_degr = QtGui.QLabel("UB_degr:", self)
		
		self.lbl_fit_out = QtGui.QLabel("fit_out:", self)
		self.lbl_fit_squ = QtGui.QLabel("fit_squ:", self)
		self.lbl_fit_slice = QtGui.QLabel("fit_slice:", self)
		self.lbl_fit_prod = QtGui.QLabel("fit_prod:", self)
		self.lbl_fit_degr = QtGui.QLabel("fit_degr:", self)
		self.lbl_eq = QtGui.QLabel("equalization:", self)
		self.lbl_pin = QtGui.QLabel("fit_pinned:", self)
		
		self.lbl_cut_off = QtGui.QLabel("Cut off at timepoint:", self)
		self.lbl_cut_off_t = QtGui.QLabel("cut_off_t:", self)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_meth = QtGui.QComboBox(self)
		self.combo_meth.addItem("TNC")
		self.combo_meth.addItem("Nelder-Mead")
		self.combo_meth.addItem("L-BFGS-B")
		self.combo_meth.addItem("SLSQP")
		self.combo_meth.addItem("brute")
		self.combo_meth.addItem("BFGS")
		self.combo_meth.addItem("CG")
		self.combo_meth.activated[str].connect(self.sel_meth)   
	
	
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		
		self.qle_opt_tol = QtGui.QLineEdit(str(self.molecule.sel_fits[0].opt_tol))
		self.qle_maxfun = QtGui.QLineEdit(str(self.molecule.sel_fits[0].maxfun))
		
		self.qle_x0_D = QtGui.QLineEdit(str(self.molecule.sel_fits[0].x0[0]))
		self.qle_x0_degr = QtGui.QLineEdit(str(self.molecule.sel_fits[0].x0[1]))
		self.qle_x0_prod = QtGui.QLineEdit(str(self.molecule.sel_fits[0].x0[2]))
		self.qle_cut_off_t = QtGui.QLineEdit(str(self.molecule.sel_fits[0].cut_off_t))
		
		if self.molecule.sel_fits[0].LB_D==None:
			self.qle_LB_D = QtGui.QLineEdit("0.0")
			self.qle_LB_D.setReadOnly(True)
		else:	
			
			self.qle_LB_D = QtGui.QLineEdit(str(self.molecule.sel_fits[0].LB_D))
			
		if self.molecule.sel_fits[0].UB_D==None:
			self.qle_UB_D = QtGui.QLineEdit("0.0")
			self.qle_UB_D.setReadOnly(True)
		else:
			self.qle_UB_D = QtGui.QLineEdit(str(self.molecule.sel_fits[0].UB_D))
			
		if self.molecule.sel_fits[0].LB_degr==None:
			self.qle_LB_degr = QtGui.QLineEdit("0.0")
			self.qle_LB_degr.setReadOnly(True)
		else:
			self.qle_LB_degr = QtGui.QLineEdit(str(self.molecule.sel_fits[0].LB_degr))
			
		if self.molecule.sel_fits[0].UB_degr==None:
			self.qle_UB_degr = QtGui.QLineEdit("0.0")
			self.qle_UB_degr.setReadOnly(True)
		else:
			self.qle_UB_degr = QtGui.QLineEdit(str(self.molecule.sel_fits[0].UB_degr))
			
		if self.molecule.sel_fits[0].LB_prod==None:
			self.qle_LB_prod = QtGui.QLineEdit("0.0")
			self.qle_LB_prod.setReadOnly(True)
		else:
			self.qle_LB_prod = QtGui.QLineEdit(str(self.molecule.sel_fits[0].LB_prod))
			
		if self.molecule.sel_fits[0].UB_prod==None:			
			self.qle_UB_prod = QtGui.QLineEdit("0.0")
			self.qle_UB_prod.setReadOnly(True)
		else:
			self.qle_UB_prod = QtGui.QLineEdit(str(self.molecule.sel_fits[0].UB_prod))
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_LB_D.setValidator(self.double_valid)
		self.qle_LB_degr.setValidator(self.double_valid)
		self.qle_LB_prod.setValidator(self.double_valid)
		self.qle_UB_D.setValidator(self.double_valid)
		self.qle_UB_degr.setValidator(self.double_valid)
		self.qle_UB_prod.setValidator(self.double_valid)
		self.qle_cut_off_t.setValidator(self.double_valid)
		
		
		self.qle_opt_tol.textChanged[str].connect(self.set_opt_tol)
		self.qle_maxfun.textChanged[str].connect(self.set_maxfun)
		
		self.qle_x0_D.textChanged[str].connect(self.set_x0_D)
		self.qle_x0_degr.textChanged[str].connect(self.set_x0_degr)
		self.qle_x0_prod.textChanged[str].connect(self.set_x0_prod)
		
		self.qle_UB_D.textChanged[str].connect(self.set_UB_D)
		self.qle_UB_degr.textChanged[str].connect(self.set_UB_degr)
		self.qle_UB_prod.textChanged[str].connect(self.set_UB_prod)
		
		self.qle_LB_D.textChanged[str].connect(self.set_LB_D)
		self.qle_LB_degr.textChanged[str].connect(self.set_LB_degr)
		self.qle_LB_prod.textChanged[str].connect(self.set_LB_prod)
		
		self.qle_cut_off_t.textChanged[str].connect(self.set_cut_off_t)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_debug_fit = QtGui.QCheckBox('', self)
		self.cb_save_track = QtGui.QCheckBox('', self)
		
		self.cb_bound_LB_D = QtGui.QCheckBox('', self)
		self.cb_bound_LB_degr = QtGui.QCheckBox('', self)
		self.cb_bound_LB_prod = QtGui.QCheckBox('', self)
		
		self.cb_bound_UB_D = QtGui.QCheckBox('', self)
		self.cb_bound_UB_degr = QtGui.QCheckBox('', self)
		self.cb_bound_UB_prod = QtGui.QCheckBox('', self)
		
		self.cb_fit_out = QtGui.QCheckBox('', self)
		self.cb_fit_squ = QtGui.QCheckBox('', self)
		self.cb_fit_slice = QtGui.QCheckBox('', self)
		
		self.cb_fit_degr = QtGui.QCheckBox('', self)
		self.cb_fit_prod = QtGui.QCheckBox('', self)
		
		self.cb_eq = QtGui.QCheckBox('', self)
		self.cb_pin = QtGui.QCheckBox('', self)
		
		self.cb_cut = QtGui.QCheckBox('', self)
		
		self.cb_debug_fit.setCheckState(2*self.molecule.sel_fits[0].debug_fit)	
		self.cb_save_track.setCheckState(2*self.molecule.sel_fits[0].save_track)	
		
		self.cb_eq.setCheckState(2*self.molecule.sel_fits[0].equ_on)	
		self.cb_pin.setCheckState(2*self.molecule.sel_fits[0].fit_pinned)	
		
		self.cb_cut.setCheckState(2*self.molecule.sel_fits[0].fit_cut_off_t)
		
		if self.molecule.sel_fits[0].LB_D==None:
			self.cb_bound_LB_D.setCheckState(0)
		else: 
			self.cb_bound_LB_D.setCheckState(2)
		
		if self.molecule.sel_fits[0].UB_D==None:
			self.cb_bound_UB_D.setCheckState(0)
		else: 
			self.cb_bound_UB_D.setCheckState(2)
		
		if self.molecule.sel_fits[0].LB_degr==None:
			self.cb_bound_LB_degr.setCheckState(0)
		else: 
			self.cb_bound_LB_degr.setCheckState(2)
		
		if self.molecule.sel_fits[0].UB_degr==None:
			self.cb_bound_UB_degr.setCheckState(0)
		else: 
			self.cb_bound_UB_degr.setCheckState(2)
		
		if self.molecule.sel_fits[0].LB_prod==None:
			self.cb_bound_LB_prod.setCheckState(0)
		else: 
			self.cb_bound_LB_prod.setCheckState(2)
		
		if self.molecule.sel_fits[0].UB_prod==None:
			self.cb_bound_UB_prod.setCheckState(0)
		else: 
			self.cb_bound_UB_prod.setCheckState(2)
		
		if self.molecule.sel_fits[0].fit_cut_off_t==1:
			self.lbl_cut_off_t.setVisible(True)
			self.qle_cut_off_t.setVisible(True)
		else: 
			self.lbl_cut_off_t.setVisible(False)
			self.qle_cut_off_t.setVisible(False)
			
		
		self.connect(self.cb_bound_LB_D, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_D)
		self.connect(self.cb_bound_LB_degr, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_degr)
		self.connect(self.cb_bound_LB_prod, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_prod)
		self.connect(self.cb_bound_UB_D, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_D)
		self.connect(self.cb_bound_UB_degr, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_degr)
		self.connect(self.cb_bound_UB_prod, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_prod)
		
		self.cb_fit_out.setCheckState(2*self.molecule.sel_fits[0].fit_out)
		self.cb_fit_squ.setCheckState(2*self.molecule.sel_fits[0].fit_squ)
		self.cb_fit_slice.setCheckState(2*self.molecule.sel_fits[0].fit_slice)
		self.cb_fit_degr.setCheckState(2*self.molecule.sel_fits[0].fit_degr)
		self.cb_fit_prod.setCheckState(2*self.molecule.sel_fits[0].fit_prod)
		
		self.connect(self.cb_debug_fit, QtCore.SIGNAL('stateChanged(int)'), self.check_debug_fit)
		self.connect(self.cb_save_track, QtCore.SIGNAL('stateChanged(int)'), self.check_save_track)
		
		
		self.connect(self.cb_fit_out, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_out)
		self.connect(self.cb_fit_squ, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_squ)
		self.connect(self.cb_fit_slice, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_slice)
		self.connect(self.cb_fit_degr, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_degr)
		self.connect(self.cb_fit_prod, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_prod)
		
		self.connect(self.cb_eq, QtCore.SIGNAL('stateChanged(int)'), self.check_equ)
		self.connect(self.cb_pin, QtCore.SIGNAL('stateChanged(int)'), self.check_pin)
		self.connect(self.cb_cut, QtCore.SIGNAL('stateChanged(int)'), self.check_cut)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_opt_parms,0,1,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_guess,0,3,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounds,0,6,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounded,0,9,1,1,Qt.AlignHCenter)
		grid.addWidget(self.lbl_fitting,0,10,1,2,Qt.AlignHCenter)
		
		

		grid.addWidget(self.lbl_opt_meth,2,1)
		grid.addWidget(self.lbl_opt_tol,3,1)
		grid.addWidget(self.lbl_maxfun,4,1)
		grid.addWidget(self.lbl_debug_fit,5,1)
		grid.addWidget(self.lbl_save_track,6,1)
		
	
		grid.addWidget(self.combo_meth,2,2)
		grid.addWidget(self.qle_opt_tol,3,2)
		grid.addWidget(self.qle_maxfun,4,2)
		grid.addWidget(self.cb_debug_fit,5,2)
		grid.addWidget(self.cb_save_track,6,2)
		
		grid.addWidget(self.lbl_x0_D,1,3)
		grid.addWidget(self.lbl_x0_degr,2,3)
		grid.addWidget(self.lbl_x0_prod,3,3)
		
		grid.addWidget(self.qle_x0_D,1,4)
		grid.addWidget(self.qle_x0_degr,2,4)
		grid.addWidget(self.qle_x0_prod,3,4)
		
		
		
		grid.addWidget(self.lbl_LB_D,1,6)
		grid.addWidget(self.lbl_LB_degr,3,6)
		grid.addWidget(self.lbl_LB_prod,5,6)
		
		grid.addWidget(self.lbl_UB_D,2,6)
		grid.addWidget(self.lbl_UB_degr,4,6)
		grid.addWidget(self.lbl_UB_prod,6,6)
		
		grid.addWidget(self.qle_LB_D,1,7)
		grid.addWidget(self.qle_LB_degr,3,7)
		grid.addWidget(self.qle_LB_prod,5,7)
		
		grid.addWidget(self.qle_UB_D,2,7)
		grid.addWidget(self.qle_UB_degr,4,7)
		grid.addWidget(self.qle_UB_prod,6,7)
	
		
		grid.addWidget(self.cb_bound_LB_D,1,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_degr,3,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_prod,5,9,Qt.AlignHCenter)
		
		grid.addWidget(self.cb_bound_UB_D,2,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_degr,4,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_prod,6,9,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_fit_out,1,10)
		grid.addWidget(self.lbl_fit_squ,2,10)
		grid.addWidget(self.lbl_fit_slice,3,10)
		grid.addWidget(self.lbl_fit_degr,4,10)
		grid.addWidget(self.lbl_fit_prod,5,10)
		grid.addWidget(self.lbl_eq,6,10)
		grid.addWidget(self.lbl_pin,7,10)
		grid.addWidget(self.lbl_cut_off,8,10)
		grid.addWidget(self.lbl_cut_off_t,9,10)
		
		grid.addWidget(self.cb_fit_out,1,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_squ,2,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_slice,3,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_degr,4,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_prod,5,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_eq,6,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_pin,7,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_cut,8,11,Qt.AlignHCenter)
		grid.addWidget(self.qle_cut_off_t,9,11,Qt.AlignHCenter)
		
		
		grid.addWidget(self.btn_done,10,11)
		
		grid.setColumnStretch(0,1)
		grid.setColumnStretch(12,1)
	
		#grid.setRowStretch(0,1)
		grid.setRowStretch(9,1)
		
		self.setLayout(grid)    
		self.setWindowTitle('Edit Fit')    
		self.show()
	
	def set_opt_tol(self,text):
		for fit in self.molecule.sel_fits:
			fit.opt_tol=float(str(text))
		
	def set_maxfun(self,text):
		for fit in self.molecule.sel_fits:
			fit.maxfun=int(str(text))
		
	def set_x0_D(self,text):
		for fit in self.molecule.sel_fits:
			fit.x0[0]=float(str(text))
		
	def set_x0_degr(self,text):
		for fit in self.molecule.sel_fits:
			fit.x0[1]=float(str(text))
			
	def set_x0_prod(self,text):
		for fit in self.molecule.sel_fits:
			fit.x0[2]=float(str(text))
			
	def set_UB_D(self,text):
		for fit in self.molecule.sel_fits:
			fit.UB_D=float(str(text))
		
	def set_UB_degr(self,text):
		for fit in self.molecule.sel_fits:
			fit.UB_degr=float(str(text))
	
	def set_UB_prod(self,text):
		for fit in self.molecule.sel_fits:
			fit.UB_prod=float(str(text))
			
	def set_LB_D(self,text):
		for fit in self.molecule.sel_fits:
			fit.LB_D=float(str(text))
		
	def set_LB_degr(self,text):
		for fit in self.molecule.sel_fits:
			fit.LB_degr=float(str(text))
	
	def set_LB_prod(self,text):
		for fit in self.molecule.sel_fits:
			fit.LB_prod=float(str(text))
	
	def set_cut_off_t(self,text):
		for fit in self.molecule.sel_fits:
			fit.cut_off_t=float(str(text))
	
	def sel_meth(self,text):
		for fit in self.molecule.sel_fits:
			fit.opt_meth=str(text)
		self.update_bounds_after_meth(str(text))
		
	def update_bounds_after_meth(self,text):
		
		
		
		if text not in ["TNC","L-BFGS-B","SLSQP","brute"]:
			
			self.temp_LB_D=fit.LB_D
			self.temp_UB_D=fit.UB_D
			self.temp_LB_degr=fit.LB_degr
			self.temp_UB_degr=fit.UB_degr
			self.temp_LB_prod=fit.LB_prod
			self.temp_UB_prod=fit.UB_prod
			
			self.check_bound_LB_D(0)
			self.check_bound_UB_D(0)
			self.check_bound_LB_degr(0)
			self.check_bound_UB_degr(0)
			self.check_bound_LB_prod(0)
			self.check_bound_UB_prod(0)
			
			self.cb_bound_LB_D.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_D.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Unchecked)
			
			
		else:
			if self.temp_LB_D==[] or text=="brute":
				
				self.check_bound_LB_D(2)
				self.check_bound_UB_D(2)
				self.check_bound_LB_degr(2)
				self.check_bound_UB_degr(2)
				self.check_bound_LB_prod(2)
				self.check_bound_UB_prod(2)
				
				self.cb_bound_LB_D.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_D.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Checked)
				
				if text=="brute":
					self.bounds_for_brute()
				
			else: 
				if self.temp_LB_D==None:
					self.check_bound_LB_D(0)
					self.cb_bound_LB_D.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_D(2)
					self.cb_bound_LB_D.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_D==None:
					self.check_bound_UB_D(0)
					self.cb_bound_UB_D.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_D(2)
					self.cb_bound_UB_D.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_LB_degr==None:
					self.check_bound_LB_degr(0)
					self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_degr(2)
					self.cb_bound_LB_degr.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_degr==None:
					self.check_bound_UB_degr(0)
					self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_degr(2)
					self.cb_bound_UB_degr.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_LB_prod==None:
					self.check_bound_LB_prod(0)
					self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_prod(2)
					self.cb_bound_LB_prod.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_prod==None:
					self.check_bound_UB_prod(0)
					self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_degr(2)
					self.cb_bound_UB_prod.setCheckState(QtCore.Qt.Checked)
		
		return
		
	def bounds_for_brute(self):
		
		for fit in self.molecule.sel_fits:
		
			if fit.LB_D>=fit.UB_D:
				print "LBs need to be larger than UBs for D, going to fix this."
				if fit.UB_D>0:
					fit.UB_D==10*fit.LB_D
				else:
					fit.UB_D=2*fit.embryo.D
				
				self.qle_UB_D.setText(str(fit.UB_D))
				
			if fit.LB_degr>=fit.UB_degr:
				print "LBs need to be larger than UBs for degr, going to fix this."
				fit.UB_degr=1.5*fit.LB_degr+0.0000001
				self.qle_UB_degr.setText(str(fit.UB_degr))
				
			if fit.LB_prod>=fit.UB_prod:
				print "LBs need to be larger than UBs for prod, going to fix this."
				fit.UB_prod=1.5*fit.LB_prod+0.0000001
				self.qle_UB_prod.setText(str(fit.UB_prod))
			
		return
	
	def check_equ(self, value):
		for fit in self.molecule.sel_fits:
			fit.equ_on=int(value/2)
	
	def check_pin(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_pinned=int(value/2)
	
	def check_cut(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_cut_off_t=int(value/2)
		self.qle_cut_off_t.setVisible(bool(int(value/2)))
		self.lbl_cut_off_t.setVisible(bool(int(value/2)))
		
	def check_fit_out(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_out=int(value/2)
		
	def check_fit_slice(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_slice=int(value/2)
			
	def check_fit_squ(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_squ=int(value/2)
				
	def check_fit_degr(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_degr=int(value/2)
	
	def check_fit_prod(self, value):
		for fit in self.molecule.sel_fits:
			fit.fit_prod=int(value/2)
	
	def check_debug_fit(self, value):
		for fit in self.molecule.sel_fits:
			fit.debug_fit=int(value/2)
		
	def check_save_track(self, value):
		for fit in self.molecule.sel_fits:
			fit.save_track=int(value/2)	
	
	def check_bound_LB_D(self,value):
		
		if value==2:
			self.qle_LB_D.setReadOnly(False)
			for fit in self.molecule.sel_fits:
				fit.LB_D=float(self.qle_LB_D.text())	
		elif value==0:	
			self.qle_LB_D.setReadOnly(True)
			for fit in self.molecule.sel_fits:
				fit.LB_D=None
			
	def check_bound_UB_D(self,value):
		if value==2:
			self.qle_UB_D.setReadOnly(False)
			for fit in self.molecule.sel_fits:
				fit.UB_D=float(self.qle_UB_D.text())	
		elif value==0:	
			self.qle_UB_D.setReadOnly(True)
			for fit in self.molecule.sel_fits:
				fit.UB_D=None		
	
	def check_bound_LB_degr(self,value):
		if value==2:
			self.qle_LB_degr.setReadOnly(False)
			for fit in self.molecule.sel_fits:
				fit.LB_degr=float(self.qle_LB_degr.text())	
		elif value==0:	
			self.qle_LB_degr.setReadOnly(True)
			for fit in self.molecule.sel_fits:
				fit.LB_degr=None
			
	def check_bound_UB_degr(self,value):
		if value==2:
			self.qle_UB_degr.setReadOnly(False)
			for fit in self.molecule.sel_fits:
				fit.UB_degr=float(self.qle_UB_degr.text())	
		elif value==0:	
			self.qle_UB_degr.setReadOnly(True)
			for fit in self.molecule.sel_fits:
				fit.UB_degr=None	
			
	def check_bound_LB_prod(self,value):
		if value==2:
			self.qle_LB_prod.setReadOnly(False)
			for fit in self.molecule.sel_fits:
				fit.LB_prod=float(self.qle_LB_prod.text())	
		elif value==0:	
			self.qle_LB_prod.setReadOnly(True)
			for fit in self.molecule.sel_fits:
				fit.LB_prod=None
			
	def check_bound_UB_prod(self,value):
		if value==2:
			self.qle_UB_prod.setReadOnly(False)
			for fit in self.molecule.sel_fits:
				fit.UB_prod=float(self.qle_UB_prod.text())	
		elif value==0:	
			self.qle_UB_prod.setReadOnly(True)
			for fit in self.molecule.sel_fits:
				fit.UB_prod=None			
			
	def done_pressed(self):
		problematic=0
		for fit in self.molecule.sel_fits:
			if fit.fit_out==0 and fit.fit_squ==0 and fit.fit_slice==0:
				print "Warning, you are fitting none of the regions in fit", fit.name, ", you need to check at least one!"
				problematic=1
		
		if problematic==0:
		
			self.done(1)
		return 

#===================================================================================================================================
#Dialog for About
#===================================================================================================================================

class about_dialog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(about_dialog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("This is PyFDAP version "+ parent.version, self)
		self.lbl_author = QtGui.QLabel("Author: Alexander Blaessle", self)
		self.lbl_website = QtGui.QLabel("Website: <a href="+parent.website+">"+parent.website+"</a>", self)
		self.connect(self.lbl_website, SIGNAL("linkActivated(QString)"), self.OpenURL) 

		self.btn_cancel=QtGui.QPushButton('Close')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.lbl_author)
		self.vbox.addWidget(self.lbl_website)		
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def OpenURL(self, URL): 
		QtGui.QDesktopServices().openUrl(QUrl(URL)) 
	
	def cancel(self):
	
		self.done(1)
	
				
#===================================================================================================================================
#Dialog for anaylze all progress
#===================================================================================================================================

class analyze_all_prog(QtGui.QDialog):
	
	def __init__(self,parent,molecule=None):
		super(analyze_all_prog,self).__init__(parent)
		
		self.molecule=molecule
		
		self.lbl_name = QtGui.QLabel("Analysing complete molecule ...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel analysis')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_analysis)	
		
		self.progressbar = QtGui.QProgressBar()
		self.progressbar.setMinimum(1)
		self.progressbar.setMaximum(len(self.molecule.embryos)*100)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.btn_cancel)
		self.vbox.addWidget(self.progressbar)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def cancel_analysis(self):
	
		self.accepted.emit()
	
	
class analyze_all_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
	prog_signal = QtCore.pyqtSignal(int,int)
	def __init__(self, molecule=None, parent=None):
		QtCore.QThread.__init__(self)
		self.molecule=molecule
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.molecule==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			i=0
			for embryo in self.molecule.embryos:
				
				embryo=pyfrp_img.analyze_dataset(embryo,signal=self.prog_signal,emb_count=i)
				i=i+1
			
				
			self.taskFinished.emit()


#===================================================================================================================================
#Dialog for simulation progress
#===================================================================================================================================

class simulation_prog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(simulation_prog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("Simulation in progress...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel simulation')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_simulation)	
		
		self.progressbar = QtGui.QProgressBar()
		self.progressbar.setMinimum(1)
		self.progressbar.setMaximum(100)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.progressbar)
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.setWindowTitle('Simulation Progress')    
		self.show()	
	
	def cancel_simulation(self):
	
		self.accepted.emit()
		
class simulation_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
	prog_signal = QtCore.pyqtSignal(int)
	
	
	def __init__(self, embryo=None, parent=None):
		QtCore.QThread.__init__(self)
		self.embryo=embryo
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.embryo==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			self.embryo=pyfrp_sim.simulate_diff_react(self.embryo,signal=self.prog_signal,emb_count=None)
			self.taskFinished.emit()
	
#===================================================================================================================================
#Dialog for fitting progress
#===================================================================================================================================

class fitting_prog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(fitting_prog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("Fitting in progress...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel fitting')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_fitting)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def cancel_fitting(self):
	
		self.accepted.emit()

class fitting_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, embryo=None, fit=None, gui=None, parent=None):
		QtCore.QThread.__init__(self)
		self.embryo=embryo
		self.fit=fit
		self.gui=gui
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.embryo==None or self.fit==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			self.embryo=pyfrp_fit.pin_conc(self.embryo)
			self.embryo=pyfrp_fit.parm_fitting(self.embryo,self.fit.fit_number,gui=self.gui)
			self.taskFinished.emit()
			
class fitting_mol_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, molecule=None, gui=None, parent=None):
		QtCore.QThread.__init__(self)
		self.molecule=molecule
		self.gui=gui
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.molecule==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			for embryo in self.molecule.embryos:
				print "Fitting all fits of embryo", embryo.name
				embryo=pyfrp_fit.pin_conc(embryo)
				for fit in embryo.fits:
					embryo=pyfrp_fit.parm_fitting(embryo,fit.fit_number,gui=self.gui)	
					print "Fitted", fit.name
			
			self.taskFinished.emit()			

#===================================================================================================================================
#Dialog for selecting fits for averaging molecule
#===================================================================================================================================

class select_fits(QtGui.QDialog):
	
	def __init__(self,molecule,single_fit,parent):
		super(select_fits,self).__init__(parent)
		
		self.molecule=molecule
		self.embr_in_right_list=[]
		self.fit_in_right_list=[]
		self.single_fit=single_fit
		
		self.btn_add=QtGui.QToolButton()
		self.btn_add.connect(self.btn_add, QtCore.SIGNAL('clicked()'), self.add_fit)
		self.btn_add.setArrowType(QtCore.Qt.RightArrow)
		
		self.btn_remove=QtGui.QToolButton()
		self.btn_remove.connect(self.btn_remove, QtCore.SIGNAL('clicked()'), self.remove_fit)
		self.btn_remove.setArrowType(QtCore.Qt.LeftArrow)
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		self.left_list=QtGui.QTreeWidget()
		self.left_list.setHeaderLabels(["Available Fits"])
		self.left_list.setColumnWidth(0,200)
		self.left_list.setColumnWidth(1,75)
		#self.left_list.currentItemChanged.connect(self.get_left_selections)
		self.left_list.itemDoubleClicked.connect(self.add_fit)
		
		
		self.right_list=QtGui.QTreeWidget()
		self.right_list.setHeaderLabels(["Selected Fits"])
		self.right_list.setColumnWidth(0,200)
		self.right_list.setColumnWidth(1,75)
		#self.right_list.currentItemChanged.connect(self.get_right_selections)
		self.right_list.itemDoubleClicked.connect(self.remove_fit)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.btn_add)
		self.vbox.addWidget(self.btn_remove)
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.left_list)
		self.hbox.addLayout(self.vbox)
		self.hbox.addWidget(self.right_list)
		
		self.vbox2 = QtGui.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btn_done)
		
		self.init_left_list()
		self.resize(400,500)
		self.setLayout(self.vbox2)
		self.show()	
	
	def init_left_list(self):
		
		for emb in self.molecule.embryos:
			
			if self.single_fit==1:
			
				#Add to left bar
				fitted=0
				if shape(emb.fits)>0:
					for fit in emb.fits:
						if shape(fit.squ_av_fitted_d)[0]>1:	
							fitted=1
				
				#Adding embryo to sidebar
				if fitted==1:
					self.curr_embr_node=QtGui.QTreeWidgetItem(self.left_list,[emb.name])
			
			else:
				self.curr_embr_node=QtGui.QTreeWidgetItem(self.left_list,[emb.name])
			
			
			#Add fits if they exists
			for fit in emb.fits:
					QtGui.QTreeWidgetItem(self.curr_embr_node,[fit.name])
			
			if self.single_fit==1:
				if fitted==1:
					self.left_list.expandItem(self.curr_embr_node)
			else:
				self.left_list.expandItem(self.curr_embr_node)

	def add_fit(self):
		self.get_left_selections()
		if self.left_list.currentItem()==None or self.left_list.currentItem().parent()==None:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			#Check if current embryo is already in right list
			if self.curr_embr in self.embr_in_right_list:
				#Check if fit is already in right list
				if self.curr_fit in self.fit_in_right_list:
					pass
				else:
					#Check if target embr already has a fit, if so, remove it
					self.get_right_ind(embryo_name=self.curr_embr.name)
					if self.single_fit==1:
						if self.curr_target_embr_node.childCount()>0:
							self.get_right_ind(fit_name=self.curr_target_embr_node.child(0).data(0,0).toString())
							self.curr_target_embr_node.takeChild(0)
							self.fit_in_right_list.remove(self.curr_target_fit)
							
					#Add to list of fits on right side	
					self.fit_in_right_list.append(self.curr_fit)
					new_node=QtGui.QTreeWidgetItem(self.curr_target_embr_node,[self.curr_fit.name])
					self.right_list.expandItem(self.curr_target_embr_node)
				
			else:
				#Add curr_embr to list
				self.embr_in_right_list.append(self.curr_embr)
				new_node=QtGui.QTreeWidgetItem(self.right_list,[self.curr_embr.name])
				
				#Add to list of fits on right side	
				self.fit_in_right_list.append(self.curr_fit)
				QtGui.QTreeWidgetItem(new_node,[self.curr_fit.name])
				self.right_list.expandItem(new_node)
		

	def remove_fit(self):
		self.get_right_selections()
		if self.right_list.currentItem()==None or self.right_list.currentItem().parent()==None:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			ind=self.curr_embr_node.indexOfChild(self.right_list.currentItem())
			self.curr_embr_node.takeChild(ind)
			self.fit_in_right_list.remove(self.curr_fit)
		
	
	def get_left_selections(self):

		#Get current embryo selection
		for emb in self.molecule.embryos:
			if self.left_list.currentItem().parent()==None: 
				#This is a embryo
				if self.left_list.currentItem().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.left_list.currentItem()
					self.curr_fit=None
					self.curr_fit_node=None
					break
			else:
				#This is a fit
				if self.left_list.currentItem().parent().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.left_list.currentItem().parent()
					for fit in self.curr_embr.fits:
						if fit.name==self.left_list.currentItem().data(0,0).toString():
							self.curr_fit=fit
							self.curr_fit_node=self.left_list.currentItem()
							break
		
	def get_right_selections(self):

		#Get current embryo selection
		for emb in self.molecule.embryos:
			if self.right_list.currentItem().parent()==None: 
				#This is a embryo
				if self.right_list.currentItem().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.right_list.currentItem()
					self.curr_fit=None
					self.curr_fit_node=None
					break
			else:
				#This is a fit
				if self.right_list.currentItem().parent().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.right_list.currentItem().parent()
					for fit in self.curr_embr.fits:
						if fit.name==self.right_list.currentItem().data(0,0).toString():
							self.curr_fit=fit
							self.curr_fit_node=self.right_list.currentItem()
							break	
						
	def get_left_ind(self,embryo_name=None,fit_name=None):	
		if embryo_name!=None:
			self.curr_target_embr_node=self.left_list.findItems(embryo_name,Qt.MatchExactly,0)
			self.curr_target_embr_node=self.curr_target_embr_node[0]
		if fit_name!=None:
			self.curr_target_fit_node=self.left_list.findItems(fit_name,Qt.MatchExactly,0)
			for fit in self.curr_embr.fits:
				if fit.name==fit_name:
					self.curr_target_fit=fit
				
	def get_right_ind(self,embryo_name=None,fit_name=None):	
		if embryo_name!=None:
			self.curr_target_embr_node=self.right_list.findItems(embryo_name,Qt.MatchExactly,0)
			self.curr_target_embr_node=self.curr_target_embr_node[0]
		if fit_name!=None:
			self.curr_target_fit_node=self.right_list.findItems(fit_name,Qt.MatchExactly,0)
			for fit in self.curr_embr.fits:
				if fit.name==fit_name:
					self.curr_target_fit=fit
					
	def done_pressed(self):
		self.molecule.sel_fits=self.fit_in_right_list
		
			
		self.done(1)	
		
