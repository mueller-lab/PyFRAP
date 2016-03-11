#! /usr/bin/python

import sys
from numpy import *
from PyQt4 import QtGui, QtCore
from pyfrp_sim_module import *
from pyfrp_img_module import *
from pyfrp_misc_module import *
import os, os.path
import matplotlib.pyplot as plt

#=====================================================================================================================================
#Defining global variables
#=====================================================================================================================================

#Data
global data_enc
global data_ft
global fn_datafolder
global data_res_nd
global conv_fact
global slice_depth_d

global man_det
global auto_det
global parms_det
global gauss_opt
global surf_opt

global offset_bleached_nd
global side_length_bleached_nd
global center_embr_nd
global radius_embr_nd

global squ_av_d
global out_av_d
global im_reg_data
global im_reg_ICs

global analyzed_flag

global nframes
global framerate
global tstart_d
global tend_d
global tvec_d

#Geometry
global geom
global cyl_radius_px
global cyl_height_px
global fish_inradius_px
global fish_outradius_px
global fish_dist_px
global fly_radius_px
global fly_zinterc_px
global frog_radius_px

#Simulation
global solver
global dt
global D
global D_data
global k_prod
global steps
global apply_data

#mesh
global mesh_gen
global cellvol_px
global cellvol_mu

#plotting
global slice_heights_px
global slice_number

global show_wire
global show_surf
global show_conc
global show_cont
global show_all

global out_wire
global out_surf
global out_conc
global out_cont
global out_all

global comp_flag

#Saving Settings
global vars_geom_names
global vars_geom_values
global vars_mesh_names
global vars_mesh_values
global vars_pde_names
global vars_pde_values
global vars_plot_names
global vars_plot_values

global save_sets
global save_data
global with_data
global done_is_pressed

#GUI flags
global geom_win_open
global mesh_win_open
global pde_win_open
global plot_win_open
global img_subwin_open
global sim_subwin_open

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Default values for global variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Data
data_enc='uint16'
data_ft='tif'
data_res_nd=512.
fn_datafolder="/home/alexander/Documents/Research/PyFRAP/Test_datasets/short_testdataset"

framerate=1

slice_depth_d=80.

slice_number=1

man_det=0
auto_det=0
parms_det=1
gauss_opt=0
surf_opt=0

offset_bleached_nd=[data_res_nd/4,data_res_nd/4]
side_length_bleached_nd=250
center_embr_nd=[data_res_nd/2,data_res_nd/2]
radius_embr_nd=300.
conv_fact=322.34/data_res_nd       
radius_data_d=conv_fact*radius_embr_nd

slice_heights_px=[-slice_depth_d/conv_fact]pyfdp_subwin.

squ_av_d=[]
out_av_d=[]
im_reg_data=[]
im_reg_ICs=[]

analyzed_flag=0

#Geometry
geom = "Cylinder"
cyl_radius_px=radius_embr_nd
cyl_height_px=4*radius_embr_nd

fish_outradius_px=(radius_embr_nd**2+slice_heights_px[0]**2)/(2*(-slice_heights_px[0]))
fish_inradius_px=fish_outradius_px*1.1
fish_dist_px=sqrt(fish_inradius_px**2-fish_outradius_px**2)


fly_radius_px=radius_embr_nd
fly_zinterc_px=0.4*radius_embr_nd

frog_radius_px=(slice_heights_px[0]**2+radius_embr_nd**2)/(2*(-slice_heights_px[0]))
			
#Simulation
D_data= 10
D=1
k_prod= 0.
steps = 4000.
dt = 1.
solver = "FiPy"
apply_data = 1

#mesh
mesh_gen = "Gmsh"
cellvol_mu = 10.
cellvol_px = cellvol_mu/conv_fact

#plotting
show_wire=0
show_surf=0
show_cont=0
show_conc=0
show_all=0

out_wire=0
out_surf=0
out_cont=0
out_conc=0
out_all=0

comp_flag=0

#Saving settings

save_sets=0
save_data=0
with_data=0
done_is_pressed=0

vars_data_names=["fn_datafolder","man_det","auto_det","parms_det","gauss_opt","surf_opt","data_enc","data_res_nd","data_ft","center_embr_nd","radius_embr_nd","offset_bleached_nd","side_length_bleached_nd","with_data"]
vars_data_values=[fn_datafolder,man_det,auto_det,parms_det,gauss_opt,surf_opt,data_enc,data_res_nd,data_ft,center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd,with_data]
vars_geom_names=["geom","cyl_radius_px","cyl_height_px"]
vars_geom_values=[geom,cyl_radius_px,cyl_height_px]
vars_mesh_names=["mesh_gen","cellvol_px"]
vars_mesh_values=[mesh_gen,cellvol_px]
vars_pde_names=["solver","steps","dt","D","k_prod"]
vars_pde_values=[solver,steps,dt,D,k_prod]
vars_plot_names=["show_wire","show_surf","show_conc","show_cont","show_all","out_wire","out_surf","out_conc","out_cont","out_all","slice_heights_px","slice_number"]
vars_plot_values=[show_wire,show_surf,show_conc,show_cont,show_all,out_wire,out_surf,out_conc,out_cont,out_all,slice_heights_px,slice_number]

#GUI Flags
geom_win_open=0
mesh_win_open=0
pde_win_open=0
plot_win_open=0
img_subwin_open=0
sim_subwin_open=0

#=====================================================================================================================================
#Main Simulation window
#=====================================================================================================================================



class pyfrp_Main(QtGui.QWidget):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		
		self.setWindowTitle('PyFRAP')

		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------

		btn_imgana=QtGui.QPushButton('Image Analysis')
		btn_sim=QtGui.QPushButton('Simulation')
		btn_fit=QtGui.QPushButton('Parameter Fitting')
		btn_stat=QtGui.QPushButton('Statistics')
		btn_save=QtGui.QPushButton('Save')
		btn_load=QtGui.QPushButton('Load')
		btn_close=QtGui.QPushButton('Close')
		
		btns=[btn_imgana,btn_sim,btn_fit,btn_stat,btn_save,btn_load,btn_close]
		
		btn_imgana.connect(btn_imgana, QtCore.SIGNAL('clicked()'), self.open_img_subwin)
		btn_sim.connect(btn_sim, QtCore.SIGNAL('clicked()'), self.open_sim_subwin)
		self.connect(btn_close, QtCore.SIGNAL('clicked()'),QtGui.qApp, QtCore.SLOT('quit()'))
		btn_save.connect(btn_save, QtCore.SIGNAL('clicked()'), self.press_save)
		btn_load.connect(btn_load, QtCore.SIGNAL('clicked()'), self.press_load)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0),  
			(1, 0),  
			(2, 0),  
			(3, 0),
			(4, 0),
			(5, 0),
			(6, 0)]
			
		#Defining Order of object types
		order=[1,1,1,1,1,1,1]
		
		#Object counters
		j=0
		k=0
		l=0
		n=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
			if order[i]==1:
				curr_obj = btns[l]
				l=l+1
				
			#Adding new obj to grid layout
			grid.addWidget(curr_obj, pos[i][0], pos[i][1])
			
		self.setLayout(grid)
	
	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------
	
	#Open Image processing subwindow
	def open_img_subwin(self):
		
		self.img_subwin=Img_Main()
		self.img_subwin.show()	
		
		global img_subwin_open
		img_subwin_open=1
		
	#Open simulation subwindow
	def open_sim_subwin(self):
		
		self.sim_subwin=Sim_Main()
		self.sim_subwin.show()	
		
		global sim_subwin_open
		sim_subwin_open=1
		
	#Load settings and data
	def press_load(self):
		
		#Loading setting file
		fn_load = QtGui.QFileDialog.getOpenFileName(self, 'Open settings file')
		settings_names,settings_vars=load_settings(fn_load)
	
		#Add to global variables
		self.add_to_globals(settings_names,settings_vars)
		
		
		#Update subwindows
		self.update_subwindows
		
		#Loading analyzed data
		if with_data==1:
			
			loaded_data=load_data_series_fast(fn_load)
			
			self.add_to_globals(["im_reg_data"],[loaded_data])
			
			global im_reg_ICs
			im_reg_ICs=im_reg_data[0]
			
			global analyzed_flag
			analyzed_flag=1
			
			
			
	#Add list of variables to globals
	def add_to_globals(self,names,varvalues):
		
		for i in range(shape(names)[0]):
			globals()[names[i]] = varvalues[i]
	
	#Updating subwindows
	def update_subwindows(self):
		
		#Updating labels and combos in subwindows if open
		if img_subwin_open==1:
			self.img_subwin.lbl_offset.setText("offset="+str(offset_bleached_nd), self)
			self.img_subwin.lbl_sidelength.setText("side length="+str(side_length_bleached_nd), self)
			self.img_subwin.lbl_center.setText("center="+str(center_embr_nd), self)
			self.img_subwin.lbl_radius.setText("radius="+str(radius_embr_nd), self)
			self.img_subwin.lbl_res.setText("resolution="+str(data_res_nd), self)
		
			self.img_subwin.update_cbs()
			
		if sim_subwin_open==1:
			if geom=="Cylinder":
				self.sim_subwin.combo_geom.setCurrentIndex(self, 0)
			elif geom=="Frog":
				self.sim_subwin.combo_geom.setCurrentIndex(self, 3)
			elif geom=="Fly":
				self.sim_subwin.combo_geom.setCurrentIndex(self, 2)
			elif geom=="Fish":
				self.sim_subwin.combo_geom.setCurrentIndex(self,1)
			
			if solver=="FiPy":
				self.sim_subwin.combo_solver.setCurrentIndex(self, 0)
			else: 
				print "Does not work yet"
				
			if mesh_gen=="Gmsh":
				self.sim_subwin.combo_mesh.setCurrentIndex(self, 0)
			else: 
				print "Does not work yet"
				
		if geom_win_open==1:
			
			if geom=="Cylinder":
				self.sim_subwin.parms_geom.lbl_radius_geom.setText("radius="+str(cyl_radius_px))
				self.sim_subwin.parms_geom.lbl_height_geom.setText("height="+str(cyl_height_px))
				self.sim_subwin.parms_geom.lbl_geometry_geom.setText("Geometry="+geom)
			
			elif geom=="Frog":          
				self.sim_subwin.lbl_radius_geom = QtGui.QLabel("radius="+str(frog_radius_px), self)
				self.sim_subwin.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
			
			elif geom=="Fly":
				self.sim_subwin.lbl_radius_geom = QtGui.QLabel("radius="+str(fly_radius_px), self)
				self.sim_subwin.lbl_zinterc_geom = QtGui.QLabel("z-intercept="+str(fly_zinterc_px), self)
				self.sim_subwin.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
			
			elif geom=="Fish":
				self.sim_subwin.lbl_outradius_geom = QtGui.QLabel("outer radius="+str(fish_outradius_px), self)
				self.sim_subwin.lbl_inradius_geom = QtGui.QLabel("inner radius="+str(fish_inradius_px), self)
				self.sim_subwin.lbl_dist_geom = QtGui.QLabel("distance between centers="+str(fish_dist_px), self)
				self.sim_subwin.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
				
		if mesh_win_open==1:
			self.sim_subwin.parms_mesh.lbl_cellvol_px.setText("Cell Volume="+str(cellvol_px))
			self.sim_subwin.parms_mesh.lbl_mesh_gen.setText("Mesh-Generator="+mesh_gen)
			
		if pde_win_open==1:
			
			self.sim_subwin.parms_pde.lbl_dt.setText("dt="+str(dt))
			self.sim_subwin.parms_pde.lbl_steps.setText("steps="+str(steps))
			self.sim_subwin.parms_pde.lbl_D.setText("D="+str(D))
			self.sim_subwin.parms_pde.lbl_k_prod.setText("k="+str(k_prod))
			
		if plot_win_open==1:
			
			#Fusing string for label
			str_slices=''
			for i in range(slice_number):
				
				if i==0:
					str_slices=str_slices+str(slice_heights_px[i])
				else:
					str_slices=str_slices+','+str(slice_heights_px[i])
			
			#Labels
			self.sim_subwin.parms_plot.lbl_nslices.setText("Number of slices="+str(slice_number))
			self.sim_subwin.parms_plot.lbl_hslices.setText("Number of slices="+str_slices)
			
			
			#Checkboxes
			if show_conc==1:
				self.sim_subwin.parms_plot.cb_show_conc.setCheckState(QtCore.Qt.Checked)
			elif show_conc==0:
				self.sim_subwin.parms_plot.cb_show_conc.setCheckState(QtCore.Qt.Unchecked)
			if show_wire==1:
				self.sim_subwin.parms_plot.cb_show_wire.setCheckState(QtCore.Qt.Checked)
			elif show_wire==0:
				self.sim_subwin.parms_plot.cb_show_wire.setCheckState(QtCore.Qt.Unchecked)
			if show_surf==1:
				self.sim_subwin.parms_plot.cb_show_surf.setCheckState(QtCore.Qt.Checked)
			elif show_surf==0:
				self.sim_subwin.parms_plot.cb_show_surf.setCheckState(QtCore.Qt.Unchecked)
			if show_cont==1:
				self.sim_subwin.parms_plot.cb_show_cont.setCheckState(QtCore.Qt.Checked)
			elif show_cont==0:
				self.sim_subwin.parms_plot.cb_show_cont.setCheckState(QtCore.Qt.Unchecked)
			if show_all==1:
				self.sim_subwin.parms_plot.cb_show_all.setCheckState(QtCore.Qt.Checked)
			elif show_all==0:
				self.sim_subwin.parms_plot.cb_show_all.setCheckState(QtCore.Qt.Unchecked)
			if out_conc==1:
				self.sim_subwin.parms_plot.cb_out_conc.setCheckState(QtCore.Qt.Checked)
			elif out_conc==0:
				self.sim_subwin.parms_plot.cb_out_conc.setCheckState(QtCore.Qt.Unchecked)
			if out_wire==1:
				self.sim_subwin.parms_plot.cb_out_wire.setCheckState(QtCore.Qt.Checked)
			elif out_wire==0:
				self.sim_subwin.parms_plot.cb_out_wire.setCheckState(QtCore.Qt.Unchecked)
			if out_surf==1:
				self.sim_subwin.parms_plot.cb_out_surf.setCheckState(QtCore.Qt.Checked)
			elif out_surf==0:
				self.sim_subwin.parms_plot.cb_out_surf.setCheckState(QtCore.Qt.Unchecked)
			if out_cont==1:
				self.sim_subwin.parms_plot.cb_out_cont.setCheckState(QtCore.Qt.Checked)
			elif out_cont==0:
				self.sim_subwin.parms_plot.cb_out_cont.setCheckState(QtCore.Qt.Unchecked)
			if out_all==1:
				self.sim_subwin.parms_plot.cb_out_all.setCheckState(QtCore.Qt.Checked)
			elif out_all==0:
				self.sim_subwin.parms_plot.cb_out_all.setCheckState(QtCore.Qt.Unchecked)
			if comp_flag==1:
				self.sim_subwin.parms_plot.cb_compare.setCheckState(QtCore.Qt.Checked)
			elif comp_flag==0:
				self.sim_subwin.parms_plot.cb_compare.setCheckState(QtCore.Qt.Unchecked)
	#Save button
	def press_save(self):
		
		self.save_set_subwin=Save_Opt_Sub()
		self.save_set_subwin.show()	
			
class Save_Opt_Sub(QtGui.QWidget):
	def __init__(self, parent=None):
		super(Save_Opt_Sub,self).__init__(parent)
		self.setWindowTitle('Save Options')
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		
		
		btns=[self.btn_done]
		
		##Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------

		self.lbl_sets = QtGui.QLabel("Save settings", self)
		self.lbl_data = QtGui.QLabel("Save analyzed data", self)
		
		lbls=[self.lbl_sets,self.lbl_data]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------

		self.cb_sets = QtGui.QCheckBox('', self)
		self.cb_data = QtGui.QCheckBox('', self)
		
		
		cbs=[self.cb_sets,self.cb_data]
		
		self.connect(self.cb_sets, QtCore.SIGNAL('stateChanged(int)'), self.check_sets)
		self.connect(self.cb_data, QtCore.SIGNAL('stateChanged(int)'), self.check_data)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0), (0, 1),  
			(1, 0), (1, 1),
			(2, 0), (2, 1)]
			
		
		#Defining Order of object types
		order=[1,3,1,3,0,2]
		
		#Object counters
		j=0
		k=0
		l=0
		m=0
		n=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
			if order[i]==1:
				curr_obj = lbls[k]
				k=k+1
				
			elif order[i]==0:
				curr_obj= QtGui.QLabel('')
				m=m+1
			
			elif order[i]==3:
				curr_obj= cbs[n]
				n=n+1
		
			elif order[i]==2:
				curr_obj = btns[l]
				l=l+1
		
		
			#Adding new obj to grid layout
			grid.addWidget(curr_obj, pos[i][0], pos[i][1])
		
		self.setLayout(grid)
		
	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------	
	
	#Save settings Checkbox
	def check_sets(self, value):
		
		global save_sets
		
		if self.cb_sets.isChecked():
			save_sets=1	
		else:
			save_sets=0
			
	#Save data Checkbox
	def check_data(self, value):
		
		global save_data
		
		if self.cb_data.isChecked():
			
			if analyzed_flag==1:
				save_data=1
			else:
				QtGui.QMessageBox.warning(self, 'Warning', "ERROR: You haven't analyzed any data yet!", QtGui.QMessageBox.Ok)
				save_data=0
				self.cb_data.setCheckState(QtCore.Qt.Unchecked)	
		else:
			save_data=0	
	
	def done_pressed(self):
			
		if save_data==1 or save_sets==1:
			#Choosing filename
			fn_save = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '.')
		
		else: 
			QtGui.QMessageBox.warning(self, 'Warning', "Nothing will be saved", QtGui.QMessageBox.Ok)
		
		if save_data==1 and save_sets==1:
			with_data=1
		else:
			with_data=0
			
			
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Settings save
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		if save_sets==1:
			vars_data_names=["fn_datafolder","man_det","auto_det","parms_det","gauss_opt","surf_opt","data_enc","data_res_nd","data_ft","center_embr_nd","radius_embr_nd","offset_bleached_nd","side_length_bleached_nd","with_data"]
			vars_data_values=[fn_datafolder,man_det,auto_det,parms_det,gauss_opt,surf_opt,data_enc,data_res_nd,data_ft,center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd,with_data]
			
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#Geometry save
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			if geom=="Cylinder":
				
				#Updating vars
				vars_geom_names=["geom","cyl_radius_px","cyl_height_px"]
				vars_geom_values=[geom,cyl_radius_px,cyl_height_px]
			
			elif geom=="Fish":
				
				#Updating vars
				vars_geom_names=["geom","fish_inradius_px","fish_outradius_px","fish_dist_px"]
				vars_geom_values=[geom,fish_inradius_px,fish_outradius_px,fish_dist_px]
				
			elif geom=="Frog":
				
				#Updating vars
				vars_geom_names=["geom","frog_radius_px"]
				vars_geom_values=[geom,frog_radius_px]
				
			elif geom=="Fly":
				
				#Updating vars
				vars_geom_names=["geom","fly_radius_px","fly_zinterc_px"]
				vars_geom_values=[geom,fly_radius_px,fly_zinterc_px]
			
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#Mesh save
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			vars_mesh_names=["mesh_gen","cellvol_px"]
			vars_mesh_values=[mesh_gen,cellvol_px]
			
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#PDE save
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			vars_pde_names=["solver","steps","dt","D","k_prod"]
			vars_pde_values=[solver,steps,dt,D,k_prod]
			
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#Plot save
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			vars_plot_names=["show_wire","show_surf","show_conc","show_cont","show_all","out_wire","out_surf","out_conc","out_cont","out_all","slice_heights_px","slice_number"]
			vars_plot_values=[show_wire,show_surf,show_conc,show_cont,show_all,out_wire,out_surf,out_conc,out_cont,out_all,slice_heights_px,slice_number]
			
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#Saving settings
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			
			save_settings(vars_data_names,vars_data_values,vars_geom_names,vars_geom_values,vars_mesh_names,vars_mesh_values,vars_pde_names,vars_pde_values,vars_plot_names,vars_plot_values,fn_save)
			
			print "Saved Settings"
		
		if save_data==1:
			
			save_data_series_fast(fn_save,im_reg_data)
			print "Saved data"
		
		self.close()
		
			
class Img_Main(QtGui.QWidget):
	def __init__(self, parent=None):
		super(Img_Main,self).__init__(parent)
		self.setWindowTitle('PyFRAP - Image Processing')
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		btn_select=QtGui.QPushButton('Select Dataset')
		btn_save=QtGui.QPushButton('Save')
		btn_load=QtGui.QPushButton('Load')
		btn_recognize=QtGui.QPushButton('Recognize')
		btn_analyze=QtGui.QPushButton('Analyze')
		btn_change_offset=QtGui.QPushButton('Change')
		btn_change_sidelength=QtGui.QPushButton('Change')
		btn_change_center=QtGui.QPushButton('Change')
		btn_change_radius=QtGui.QPushButton('Change')
		btn_change_res=QtGui.QPushButton('Change')
		btn_change_depth=QtGui.QPushButton('Change')
		btn_change_conv=QtGui.QPushButton('Change')
		btn_change_framerate=QtGui.QPushButton('Change')
		btns=[btn_select,btn_change_res,btn_change_conv,btn_change_offset,btn_change_sidelength,btn_change_center,btn_change_radius,btn_change_depth,btn_change_framerate,btn_recognize,btn_analyze]
		
		##Button Actions
		btn_select.connect(btn_select, QtCore.SIGNAL('clicked()'), self.select_dataset)
		
		btn_change_offset.connect(btn_change_offset, QtCore.SIGNAL('clicked()'), self.change_offset)
		btn_change_sidelength.connect(btn_change_sidelength, QtCore.SIGNAL('clicked()'), self.change_sidelength)
		btn_change_center.connect(btn_change_center, QtCore.SIGNAL('clicked()'), self.change_center)
		btn_change_radius.connect(btn_change_radius, QtCore.SIGNAL('clicked()'), self.change_radius)
		btn_change_res.connect(btn_change_res, QtCore.SIGNAL('clicked()'), self.change_res)
		btn_recognize.connect(btn_recognize, QtCore.SIGNAL('clicked()'), self.recogn_embr_bleached)
		btn_analyze.connect(btn_analyze, QtCore.SIGNAL('clicked()'), self.analyze_dataset)
		btn_change_conv.connect(btn_change_conv, QtCore.SIGNAL('clicked()'), self.change_conv)
		btn_change_depth.connect(btn_change_depth, QtCore.SIGNAL('clicked()'), self.change_depth)
		btn_change_framerate.connect(btn_change_framerate, QtCore.SIGNAL('clicked()'), self.change_framerate)
		#btn_save.connect(btn_mesh, QtCore.SIGNAL('clicked()'), self.save_img_settings)
		#btn_load.connect(btn_pde, QtCore.SIGNAL('clicked()'), self.load_img_settings)

		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
	
		
		self.lbl_mode = QtGui.QLabel("Geometry recognition mode", self)
		self.lbl_manual = QtGui.QLabel("Manual", self)
		self.lbl_auto = QtGui.QLabel("Automatic", self)
		self.lbl_parms = QtGui.QLabel("Static parameters", self)
		self.lbl_showparms = QtGui.QLabel("Parameters", self)
		self.lbl_offset = QtGui.QLabel("offset="+str(offset_bleached_nd[0]*conv_fact)+','+str(offset_bleached_nd[1]*conv_fact)+"um", self)
		self.lbl_sidelength = QtGui.QLabel("side length="+str(side_length_bleached_nd*conv_fact)+"um", self)
		self.lbl_center = QtGui.QLabel("center="+str(center_embr_nd[0]*conv_fact)+","+str(center_embr_nd[1]*conv_fact)+ "um", self)
		self.lbl_radius = QtGui.QLabel("radius="+str(radius_embr_nd*conv_fact)+ "um", self)
		self.lbl_res = QtGui.QLabel("resolution="+str(data_res_nd*conv_fact)+"um", self)
		self.lbl_surf = QtGui.QLabel("Surface Plot", self)
		self.lbl_gauss = QtGui.QLabel("Gaussian", self)
		self.lbl_conv = QtGui.QLabel("Conversion Factor="+str(conv_fact)+"px/um", self)
		self.lbl_depth = QtGui.QLabel("Slice Depth="+str(slice_depth_d), self)
		self.lbl_nframes = QtGui.QLabel("Frames=", self)
		self.lbl_framerate = QtGui.QLabel("Framerate="+str(framerate), self)
		
		self.lbl_offset_px = QtGui.QLabel(str(offset_bleached_nd[0])+','+str(offset_bleached_nd[1])+"px", self)
		self.lbl_sidelength_px = QtGui.QLabel(str(side_length_bleached_nd)+"px", self)
		self.lbl_center_px = QtGui.QLabel(str(center_embr_nd[0])+","+str(center_embr_nd[1])+ "px", self)
		self.lbl_radius_px = QtGui.QLabel(str(radius_embr_nd)+ "px", self)
		self.lbl_res_px = QtGui.QLabel(str(data_res_nd)+"px", self)
		self.lbl_depth_px = QtGui.QLabel(str(slice_depth_d/conv_fact)+"px", self)
		
		
		lbls=[self.lbl_mode,self.lbl_manual,self.lbl_auto,self.lbl_parms,self.lbl_showparms,self.lbl_res,self.lbl_res_px, self.lbl_conv,self.lbl_offset,self.lbl_offset_px,self.lbl_sidelength,self.lbl_sidelength_px,self.lbl_center,self.lbl_center_px,self.lbl_radius,self.lbl_radius_px,self.lbl_depth,self.lbl_depth_px,self.lbl_nframes,self.lbl_framerate,self.lbl_gauss,self.lbl_surf]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------

		self.cb_manual = QtGui.QCheckBox('', self)
		self.cb_auto = QtGui.QCheckBox('', self)
		self.cb_parms = QtGui.QCheckBox('', self)
		self.cb_gauss = QtGui.QCheckBox('', self)
		self.cb_surf = QtGui.QCheckBox('', self)
		
		cbs=[self.cb_manual,self.cb_auto,self.cb_parms,self.cb_gauss,self.cb_surf]
		
		self.connect(self.cb_manual, QtCore.SIGNAL('stateChanged(int)'), self.check_man_detect)
		self.connect(self.cb_auto, QtCore.SIGNAL('stateChanged(int)'), self.check_auto_detect)
		self.connect(self.cb_parms, QtCore.SIGNAL('stateChanged(int)'), self.check_parms_detect)
		self.connect(self.cb_gauss, QtCore.SIGNAL('stateChanged(int)'), self.check_gaussian)
		self.connect(self.cb_surf, QtCore.SIGNAL('stateChanged(int)'), self.check_surf)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		#filetype
		combo_ft = QtGui.QComboBox(self)
		combo_ft.addItem("TIFF,16bit")
		combo_ft.addItem("TIFF,8bit")
		
		
		combos=[combo_ft]
		
		combo_ft.activated[str].connect(self.sel_ft)   
			
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid_plot = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0), (0, 1), (0,2),  
			(1, 0), (1, 1), (1,2), 
			(2, 0), (2, 1), (2,2),
			(3, 0), (3, 1), (3,2),
			(4, 0), (4, 1), (4,2),
			(5, 0), (5, 1), (5,2),
			(6, 0), (6, 1), (6,2),
			(7, 0), (7, 1), (7,2),
			(8, 0), (8, 1), (8,2),
			(9, 0), (9, 1), (9,2),
			(10, 0), (10, 1), (10,2),
			(11, 0), (11, 1), (11,2),
			(12, 0), (12, 1), (12,2),
			(13, 0), (13, 1), (13,2),
			(14, 0), (14, 1), (14,2),
			(15, 0), (15, 1), (15,2),
			(16, 0), (16, 1), (16,2),
			(17, 0), (17, 1), (17,2)]
			
		
		#Defining Order of object types
		order=[2,4,0,1,0,0,1,3,0,1,3,0,1,3,0,1,0,0,1,1,2,1,0,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,0,0,1,0,2,1,3,0,1,3,0,2,2,0]
		
		#Object counters
		j=0
		k=0
		l=0
		m=0
		n=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
			if order[i]==1:
				curr_obj = lbls[k]
				k=k+1
	
				
			elif order[i]==2:
				curr_obj = btns[l]
				l=l+1
			
			elif order[i]==0:
				curr_obj= QtGui.QLabel('')
				m=m+1
			
			elif order[i]==3:
				curr_obj= cbs[n]
				n=n+1
		
			elif order[i]==4:
				curr_obj= combos[j]
				j=j+1
		
		
			#Adding new obj to grid layout
			grid_plot.addWidget(curr_obj, pos[i][0], pos[i][1])
		
		self.setLayout(grid_plot)

	
	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------

	#Combo for datatype
	def sel_ft(self,text):
		global data_enc
		global data_ft
		if text=="TIFF,8bit":
		
			data_ft="tif"
			data_enc="uint8"
			
		if text=="TIFF,16bit":
		
			data_ft="tif"
			data_enc="uint16"
	
	#Select data folder
	def select_dataset(self):
		
		global fn_datafolder
		fn_datafolder = QtGui.QFileDialog.getExistingDirectory(self, 'Open Data folder')
	        self.get_file_list(fn_datafolder)
	        
	##Select data folder
	def get_file_list(self,fn_datafolder):
		
		#Getting files in datafolder
		file_list=os.listdir(fn_datafolder)
		
		file_list_new=[]
		
		#Going through all 
		for i in range(shape(file_list)[0]):
			#Check if its the right picture type
			if data_ft in file_list[i]:
				file_list_new.append(file_list[i])
		
		#print shape(file_list_new)[0], " files to be loaded"
		file_list_new.sort()
	
		global nframes
		nframes= shape(file_list_new)[0]
	        self.lbl_nframes.setText("Frames="+str(nframes))
		
		global tstart_d
		global tend_d
		
		tstart_d=0
		tend_d=nframes*(framerate-1)
		
		global tvec_d
		tvec_d=linspace(tstart_d,tend_d,nframes)
		
		return file_list_new
	
	#Manual Checkbox
	def check_man_detect(self, value):
		
		global man_det
		global auto_det
		global parms_det
		
		if self.cb_manual.isChecked():
			man_det=1
			auto_det=0
			parms_det=0
			self.update_cbs()
					
		else:
			man_det=0	
			self.update_cbs()
	
	#Change button offset
	def change_offset(self):
		
		global man_det
		global auto_det
		global parms_det
		
		#Change to manual mode
		man_det=1
		auto_det=0
		parms_det=0
	
		self.update_cbs()
		
		#Dialogs
		offset_bleached_nd_x, ok = QtGui.QInputDialog.getDouble(self, 'Bleached offset x', 'offset_x=')
		offset_bleached_nd_y, ok = QtGui.QInputDialog.getDouble(self, 'Bleached offset y', 'offset_y=')
		
		global offset_bleached_nd
		offset_bleached_nd=[offset_bleached_nd_x,offset_bleached_nd_y]
		
		#Update labels
		self.lbl_offset.setText("offset="+str(offset_bleached_nd*conv_fact),"mu")
	        self.lbl_offset_px.setText(str(offset_bleached_nd)+"px")
	#Change button sidelength
	def change_sidelength(self):
		
		global man_det
		global auto_det
		global parms_det
		
		#Change to manual mode
		man_det=1
		auto_det=0
		parms_det=0
	
		self.update_cbs()
		
		#Dialogs
		global side_length_bleached_nd
		side_length_bleached_nd, ok = QtGui.QInputDialog.getDouble(self, 'Side length bleached area', 'l_sq=')
		
		#Update labels
		self.lbl_sidelength_px.setText(str(side_length_bleached_nd+"px"))
		self.lbl_sidelength.setText("side length="+str(side_length_bleached_nd*conv_fact)+"um")
		
		
	
	#Change button center
	def change_center(self):
		
		global man_det
		global auto_det
		global parms_det
		
		#Change to manual mode
		man_det=1
		auto_det=0
		parms_det=0
	
		self.update_cbs()
		
		#Dialogs
		center_embr_nd_x, ok = QtGui.QInputDialog.getDouble(self, 'Center embryo x', 'center_x=')
		center_embr_nd_y, ok = QtGui.QInputDialog.getDouble(self, 'Center embryo y', 'center_y=')
		
		global center_embr_nd
		center_embr_nd=[center_embr_nd_x,center_embr_nd_y]
		
		#Update labels
		self.lbl_center.setText("center="+str(center_embr_nd[0]*conv_fact)+','+str(center_embr_nd[1]*conv_fact)+"um")
		self.lbl_center_px.setText(str(center_embr_nd[0])+','+str(center_embr_nd[1])+"px")
	
	#Change button radius
	def change_radius(self):
		
		global man_det
		global auto_det
		global parms_det
		
		#Change to manual mode
		man_det=1
		auto_det=0
		parms_det=0
	
		self.update_cbs()
		
		#Dialogs
		global radius_embr_nd
		radius_embr_nd, ok= QtGui.QInputDialog.getDouble(self, 'embryo radius', 'r=')
		
		#Update labels
		self.lbl_radius.setText("radius="+str(radius_embr_nd*conv_fact)+ "um")
		self.lbl_radius_px.setText(str(radius_embr_nd*conv_fact)+ "px")
		
		
	#Change button radius
	def change_conv(self):
		
		#Dialogs
		global conv_fact
		conv_fact, ok= QtGui.QInputDialog.getDouble(self, 'conversion factor', 'conv_fact=')
		
		#Update labels
		self.lbl_conv.setText("conversion factor="+str(conv_fact))
	
	#Change button radius
	def change_framerate(self):
		
		#Dialogs
		global framerate
		framerate, ok= QtGui.QInputDialog.getInt(self, 'Framerate', 'framerate=')
		
		#Update labels
		self.lbl_framerate.setText("framerate="+str(framerate))
	
	#Change button depth
	def change_depth(self):
		
		#Dialogs
		global slice_depth_d
		slice_depth_d, ok= QtGui.QInputDialog.getDouble(self, 'slice_depth_d', 'z=')
		
		#Update labels
		self.lbl_depth.setText("slice depth="+str(slice_depth_d))	
		
	#Change button radius
	def change_res(self):
		
		global man_det
		global auto_det
		global parms_det
		
		#Change to manual mode
		man_det=0
		auto_det=0
		parms_det=1
	
		self.update_cbs()
		
		#Dialogs
		global data_res_nd
		data_res_nd, ok= QtGui.QInputDialog.getDouble(self, 'Set data resolution', 'res=')
		
		#Update labels
		self.lbl_res.setText("res="+str(data_res_nd*conv_fact))
	        self.lbl_res_px.setText(str(data_res_nd))
	        
	#auto Checkbox
	def check_auto_detect(self, value):
		
		global man_det
		global auto_det
		global parms_det
		
		
		if self.cb_auto.isChecked():
			man_det=0
			auto_det=1
			parms_det=0
			self.update_cbs()
		else:
			auto_det=0
			self.update_cbs()
	
	#Gaussian Checkbox
	def check_gaussian(self, value):
		
		global gauss_opt
		
		if self.cb_gauss.isChecked():
			gauss_opt=1
		else:
			gauss_opt=0
	
	#Surface Checkbox
	def check_surf(self, value):
		
		global surf_opt
		
		if self.cb_surf.isChecked():
			surf_opt=1
		else:
			surf_opt=0
	
	#Parms Checkbox
	def check_parms_detect(self, value):
		
		global man_det
		global auto_det
		global parms_det
		
		if self.cb_parms.isChecked():
			man_det=0
			auto_det=0
			parms_det=1
			
			
			self.update_cbs()
		else:
			parms_det=0	
			self.update_cbs()
	
	#Update Checkboxes
	def update_cbs(self):
		
		#manual
		if man_det==1:
			self.cb_manual.setCheckState(QtCore.Qt.Checked)
		elif man_det==0:
			self.cb_manual.setCheckState(QtCore.Qt.Unchecked)	
		#manual
		if auto_det==1:
			self.cb_auto.setCheckState(QtCore.Qt.Checked)
		elif auto_det==0:
			self.cb_auto.setCheckState(QtCore.Qt.Unchecked)	
		#manual
		if parms_det==1:
			self.cb_parms.setCheckState(QtCore.Qt.Checked)
		elif parms_det==0:
			self.cb_parms.setCheckState(QtCore.Qt.Unchecked)
			
		#Gaussian
		if gauss_opt==1:
			self.cb_gauss.setCheckState(QtCore.Qt.Checked)
		elif gauss_opt==0:
			self.cb_gauss.setCheckState(QtCore.Qt.Unchecked)
			
		#Surf
		if surf_opt==1:
			self.cb_surf.setCheckState(QtCore.Qt.Checked)
		elif surf_opt==0:
			self.cb_surf.setCheckState(QtCore.Qt.Unchecked)
			
	
	#Recognize embryo radius and bleached area
	def recogn_embr_bleached(self):
		
		#Finding stuff
		global offset_bleached_nd
		global side_length_bleached_nd
		global center_embr_nd
		global radius_embr_nd
		global data_res_nd
		
		#Parms vec
		parms_for_det=[center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd]
		
		#File list
		file_list=self.get_file_list(fn_datafolder)
		
		center_embr_nd, radius_embr_nd, offset_bleached_nd, side_length_bleached_nd, data_res_nd =find_img_bnds(fn_datafolder,file_list,data_enc,man_det,auto_det,parms_det,gauss_opt,1,surf_opt,parms_for_det)
		
		#Updating labels
		self.lbl_offset.setText("offset="+str(offset_bleached_nd))
		self.lbl_radius.setText("radius="+str(radius_embr_nd))
		self.lbl_center.setText("center="+str(center_embr_nd))
		self.lbl_sidelength.setText("side length="+str(side_length_bleached_nd))
		self.lbl_res.setText("res="+str(data_res_nd))
		self.lbl_conv.setText("Conversion Factor="+str(conv_fact))
		self.lbl_depth.setText("Conversion Factor="+str(conv_fact))
		
	#Analyze data set
	def analyze_dataset(self):
		
		#File list
		file_list=self.get_file_list(fn_datafolder)
		
		global squ_av_d
		global out_av_d
		global im_reg_data
		global im_reg_ICs
		
		squ_av_d, out_av_d, im_reg_data, im_reg_ICs = analyze_dataset(fn_datafolder,file_list,data_enc,offset_bleached_nd,side_length_bleached_nd,radius_embr_nd,center_embr_nd)
		
		global analyzed_flag
		analyzed_flag=1
		
class Sim_Main(QtGui.QWidget):
	def __init__(self, parent=None):
		super(Sim_Main,self).__init__(parent)
		self.setWindowTitle('PyFRAP - Simulation')
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------

		btn_start=QtGui.QPushButton('Start')
		btn_geom=QtGui.QPushButton('Geometry Parms')
		btn_mesh=QtGui.QPushButton('Mesh Parms')
		btn_pde=QtGui.QPushButton('PDE Parms')
		btn_plot=QtGui.QPushButton('Plot Options')
		btn_clfig=QtGui.QPushButton('Close all figures')
		
		btns=[btn_geom,btn_mesh,btn_pde,btn_plot,btn_start,btn_clfig]
		
		#Button Actions
		btn_geom.connect(btn_geom, QtCore.SIGNAL('clicked()'), self.set_geom)
		btn_mesh.connect(btn_mesh, QtCore.SIGNAL('clicked()'), self.set_mesh)
		btn_pde.connect(btn_pde, QtCore.SIGNAL('clicked()'), self.set_pde)
		btn_plot.connect(btn_plot, QtCore.SIGNAL('clicked()'), self.set_plot)
		btn_start.connect(btn_start, QtCore.SIGNAL('clicked()'), self.start_sim)
		btn_clfig.connect(btn_clfig, QtCore.SIGNAL('clicked()'), self.cl_fig)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------

		lbl_geom = QtGui.QLabel("Geometry", self)
		#lbl_mesh = QtGui.QLabel("Mesh", self)
		#lbl_solver = QtGui.QLabel("Solver", self)		
		
		lbls=[lbl_geom]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		#Geometry
		combo_geom = QtGui.QComboBox(self)
		combo_geom.addItem("Cylinder")
		combo_geom.addItem("Fish")
		combo_geom.addItem("Fly")
		combo_geom.addItem("Frog")
		
		#Mesh
		#combo_mesh = QtGui.QComboBox(self)
		#combo_mesh.addItem("Gmsh")
		#combo_mesh.addItem("Other")
		
		#Solver
		#combo_solver = QtGui.QComboBox(self)
		#combo_solver.addItem("FiPy")		
		#combo_solver.addItem("Other")
		
		#Combo Actions
		combo_geom.activated[str].connect(self.sel_geom)        
		#combo_mesh.activated[str].connect(self.sel_mesh)        
		#combo_solver.activated[str].connect(self.sel_solver)        
			
		#combos=[combo_geom,combo_mesh,combo_solver]
		combos=[combo_geom]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0), (0, 1), 
			(1, 0), (1, 1), 
			(2, 0), (2, 1), 
			(3, 0), (3, 1), 
			(4, 0), (4, 1),
			(5, 0), (5, 1)]
		
		#Defining Order of object types
		order=[1,2,0,2,3,2,3,2,3,2,3,2]
		
		#Object counters
		j=0
		k=0
		l=0
		n=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
			if order[i]==0:
				curr_obj = combos[j]
				j=j+1
				
			elif order[i]==1:
				curr_obj = lbls[k]
				k=k+1
				
			elif order[i]==2:
				curr_obj = btns[l]
				l=l+1
			
			elif order[i]==3:
				curr_obj = QtGui.QLabel('')
				n=n+1
				
			#Adding new obj to grid layout
			grid.addWidget(curr_obj, pos[i][0], pos[i][1])
			
			
			
		self.setLayout(grid)

	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------
	
	def cl_fig(self):
		plt.close('all')
		
	#Selects solver	
	def sel_solver(self, text):
		global solver	
		solver=text
		return solver
	
	#Selects geometry
	def sel_geom(self, text):	
		
		global geom	
		geom=text
		
		global slice_heights_px
		
		if geom=="Cylinder":
			cyl_radius_px=radius_embr_nd
			cyl_height_px=4*radius_embr_nd
			slice_heights_px=[-slice_depth_d/conv_fact]
			
		elif geom=="Fish":
			
			slice_heights_px=[-slice_depth_d/conv_fact]
			fish_outradius_px=(radius_embr_nd**2+slice_heights_px[0]**2)/(2*(-slice_heights_px[0]))
			fish_inradius_px=fish_outradius_px*1.1
			fish_dist_px=sqrt(fish_inradius_px**2-fish_outradius_px**2)

			
		elif geom=="Frog":
			
			slice_heights_px=[-slice_depth_d/conv_fact]
			frog_radius_px=(slice_heights_px[0]**2+radius_embr_nd**2)/(2*slice_heights_px[0])
			
		elif geom=="Fly":
			
			slice_heights_px=[-slice_depth_d/conv_fact]
			
			print "Fly Geometry is not working automatically yet."
		
		return geom 
	
	#Selects mesh generator
	def sel_mesh(self, text):
		
		global mesh_gen
		mesh_gen=text
		return mesh_gen
	
	#Opens subwindow for geometry parameters
	def set_geom(self):
		
		#Flag that window is open
		global geom_win_open
		geom_win_open=1
		
		self.parms_geom=win_parms_geom()
		self.parms_geom.show()
		
	#Opens subwindow for mesh parameters
	def set_mesh(self):
		
		#Flag that window is open
		global mesh_win_open
		mesh_win_open=1
		
		self.parms_mesh=win_parms_mesh()
		self.parms_mesh.show()	
	
	#Opens subwindow for pde parameters
	def set_pde(self):
		
		#Flag that window is open
		global pde_win_open
		pde_win_open=1
		
		self.parms_pde=win_parms_pde()
		self.parms_pde.show()	
		
	#Opens subwindow for pde parameters
	def set_plot(self):
		
		#Flag that window is open
		global plot_win_open
		plot_win_open=1
		
		self.parms_plot=win_parms_plot()
		self.parms_plot.show()	
		
	#Start simulation
	def start_sim(self):
		
		if geom=="Cylinder":
			geom_parms=[cyl_radius_px,cyl_height_px,center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd]
		elif geom=="Fish":
			geom_parms=[fish_outradius_px,fish_inradius_px,fish_dist_px,center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd]
		elif geom=="Frog":
			geom_parms=[frog_radius_px,center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd]
		elif geom=="Fly":
			geom_parms=[fly_radius_px,fly_zinterc_px,center_embr_nd,radius_embr_nd,offset_bleached_nd,side_length_bleached_nd]
		
		global sim_res_phi_reg
		global sim_res_squ_av
		global sim_res_out_av
		global sim_res_squ_small_av
		global sim_res_out_small_av
		
		
		sim_res_phi_reg, sim_res_squ_av, sim_res_out_av, sim_res_squ_small_av, sim_res_out_small_av, phi_IC_rad=D_reg(show_wire,show_surf,show_conc,show_cont,show_all,D,k_prod,dt,steps,geom,geom_parms,cellvol_px,slice_heights_px,out_wire,out_surf,out_conc,out_cont,out_all,apply_data,im_reg_ICs,0,[],0,[],0)

		
			
#=====================================================================================================================================
#Window for Geometry parameter setting 
#=====================================================================================================================================


class win_parms_geom(QtGui.QWidget):
	def __init__(self, parent=None):
		super(win_parms_geom,self).__init__(parent)
	
		self.setWindowTitle('Geometry Parameter')
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Cylinder
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		if geom=="Cylinder":
	
			#-------------------------------------------------------------------------------------------------------------------
			#Buttons
			#-------------------------------------------------------------------------------------------------------------------
			
			btn_change_radius_geom=QtGui.QPushButton('Change',self)
			btn_change_height_geom=QtGui.QPushButton('Change',self)

			btns_geom=[btn_change_radius_geom,btn_change_height_geom]
			
			#Button Actions
			btn_change_radius_geom.connect(btn_change_radius_geom, QtCore.SIGNAL('clicked()'), self.change_cyl_radius_px)
			btn_change_height_geom.connect(btn_change_height_geom, QtCore.SIGNAL('clicked()'), self.change_cyl_height_px)
				
			#-------------------------------------------------------------------------------------------------------------------
			#Labels
			#-------------------------------------------------------------------------------------------------------------------

			self.lbl_radius_geom = QtGui.QLabel("radius="+str(cyl_radius_px), self)
			self.lbl_height_geom = QtGui.QLabel("height="+str(cyl_height_px), self)
			self.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
			
			lbls_geom=[self.lbl_radius_geom,self.lbl_height_geom,self.lbl_geometry_geom]
			
			#-------------------------------------------------------------------------------------------------------------------
			#Grid
			#-------------------------------------------------------------------------------------------------------------------
			
			#Creating Grid
			grid_geom = QtGui.QGridLayout()
			
			#Defining Grid Positions
			pos = [(0, 0), (0, 1), 
				(1, 0), (1, 1), 
				(2, 0), (2, 1)] 
				
			#Defining Order of object types
			order=[1,2,1,2,1,3]
			
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Fish
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		elif geom=="Fish":
		
	
			#-------------------------------------------------------------------------------------------------------------------
			#Buttons
			#-------------------------------------------------------------------------------------------------------------------
			
			btn_change_outradius_geom=QtGui.QPushButton('Change',self)
			btn_change_inradius_geom=QtGui.QPushButton('Change',self)
			btn_change_dist_geom=QtGui.QPushButton('Change',self)
			
			btns_geom=[btn_change_outradius_geom,btn_change_inradius_geom,btn_change_dist_geom]
			
			#Button Actions
			btn_change_outradius_geom.connect(btn_change_outradius_geom, QtCore.SIGNAL('clicked()'), self.change_fish_outradius_px)
			btn_change_inradius_geom.connect(btn_change_inradius_geom, QtCore.SIGNAL('clicked()'), self.change_fish_inradius_px)
			btn_change_dist_geom.connect(btn_change_dist_geom, QtCore.SIGNAL('clicked()'), self.change_fish_dist_px)
		
				
			#-------------------------------------------------------------------------------------------------------------------
			#Labels
			#-------------------------------------------------------------------------------------------------------------------

			self.lbl_outradius_geom = QtGui.QLabel("outer radius="+str(fish_outradius_px), self)
			self.lbl_inradius_geom = QtGui.QLabel("inner radius="+str(fish_inradius_px), self)
			self.lbl_dist_geom = QtGui.QLabel("distance between centers="+str(fish_dist_px), self)
			self.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
			
			lbls_geom=[self.lbl_outradius_geom,self.lbl_inradius_geom,self.lbl_dist_geom,self.lbl_geometry_geom]
			
			#-------------------------------------------------------------------------------------------------------------------
			#Grid
			#-------------------------------------------------------------------------------------------------------------------
			
			#Creating Grid
			grid_geom = QtGui.QGridLayout()
			
			#Defining Grid Positions
			pos = [(0, 0), (0, 1), 
				(1, 0), (1, 1), 
				(2, 0), (2, 1),
				(3, 0), (3, 1)] 
				
			#Defining Order of object types
			order=[1,2,1,2,1,2,1,3]
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Fly
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		elif geom=="Fly":
		
	
			#-------------------------------------------------------------------------------------------------------------------
			#Buttons
			#-------------------------------------------------------------------------------------------------------------------
			
			btn_change_radius_geom=QtGui.QPushButton('Change',self)
			btn_change_zinterc_geom=QtGui.QPushButton('Change',self)
			
			btns_geom=[btn_change_radius_geom,btn_change_zinterc_geom]
			
			#Button Actions
			btn_change_radius_geom.connect(btn_change_radius_geom, QtCore.SIGNAL('clicked()'), self.change_fly_radius_px)
			btn_change_zinterc_geom.connect(btn_change_zinterc_geom, QtCore.SIGNAL('clicked()'), self.change_fly_zinterc_px)
				
			#-------------------------------------------------------------------------------------------------------------------
			#Labels
			#-------------------------------------------------------------------------------------------------------------------

			self.lbl_radius_geom = QtGui.QLabel("radius="+str(fly_radius_px), self)
			self.lbl_zinterc_geom = QtGui.QLabel("z-intercept="+str(fly_zinterc_px), self)
			self.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
			
			lbls_geom=[self.lbl_radius_geom,self.lbl_zinterc_geom,self.lbl_geometry_geom]
			
			#-------------------------------------------------------------------------------------------------------------------
			#Grid
			#-------------------------------------------------------------------------------------------------------------------
			
			#Creating Grid
			grid_geom = QtGui.QGridLayout()
			
			#Defining Grid Positions
			pos = [(0, 0), (0, 1), 
				(1, 0), (1, 1), 
				(2, 0), (2, 1)]
				
			#Defining Order of object types
			order=[1,2,1,2,1,3]
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Frog
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		elif geom=="Frog":
		
	
			#-------------------------------------------------------------------------------------------------------------------
			#Buttons
			#-------------------------------------------------------------------------------------------------------------------
			
			btn_change_radius_geom=QtGui.QPushButton('Change',self)
			
			btns_geom=[btn_change_radius_geom]
			
			#Button Actions
			btn_change_radius_geom.connect(btn_change_radius_geom, QtCore.SIGNAL('clicked()'), self.change_frog_radius_px)
		
				
			#-------------------------------------------------------------------------------------------------------------------
			#Labels
			#-------------------------------------------------------------------------------------------------------------------

			self.lbl_radius_geom = QtGui.QLabel("radius="+str(frog_radius_px), self)
			self.lbl_geometry_geom = QtGui.QLabel("Geometry="+geom, self)
			
			lbls_geom=[self.lbl_radius_geom,self.lbl_geometry_geom]
			
			#-------------------------------------------------------------------------------------------------------------------
			#Grid
			#-------------------------------------------------------------------------------------------------------------------
			
			#Creating Grid
			grid_geom = QtGui.QGridLayout()
			
			#Defining Grid Positions
			pos = [(0, 0), (0, 1), 
				(1, 0), (1, 1)] 
			 
				
			#Defining Order of object types
			order=[1,2,1,3]
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Building Grid
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		#Object counters
		j=0
		k=0
		l=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
				
			if order[i]==1:
				curr_obj = lbls_geom[k]
				k=k+1
				
			elif order[i]==2:
				curr_obj = btns_geom[l]
				l=l+1
			
			elif order[i]==3:
				curr_obj = QtGui.QLabel('')
				
			#Adding new obj to grid layout
			grid_geom.addWidget(curr_obj, pos[i][0], pos[i][1])
		
		self.setLayout(grid_geom)
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Action Functions
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
	#Change radius of cylinder
	def change_cyl_radius_px(self):
		new_r, ok = QtGui.QInputDialog.getDouble(self, 'Change radius', 'r=')
		
		global cyl_radius_px
		cyl_radius_px=new_r
		
		if ok:
			self.lbl_radius_geom.setText("radius="+str(new_r))
		return new_r
	
	#Change height of cylinder
	def change_cyl_height_px(self):
		new_h, ok = QtGui.QInputDialog.getDouble(self, 'Change height', 'h=')
		
		global cyl_height_px
		cyl_height_px=new_h
		
		if ok:
			self.lbl_height_geom.setText("height"+str(new_h))
		return new_h
	
	#Change outer radius of fish
	def change_fish_outradius_px(self):
		new_r, ok = QtGui.QInputDialog.getDouble(self, 'Change outer radius', 'r=')
		
		global fish_outradius_px
		fish_outradius_px=new_r
		
		if ok:
			self.lbl_outradius_geom.setText("outer radius="+str(new_r))
		return new_r
	
	#Change inradius radius of fish
	def change_fish_inradius_px(self):
		new_r, ok = QtGui.QInputDialog.getDouble(self, 'Change inner radius', 'r=')
		
		global fish_inradius_px
		fish_inradius_px=new_r
		
		if ok:
			self.lbl_inradius_geom.setText("inner radius="+str(new_r))
		return new_r
	
	#Change center_dist
	def change_fish_dist_px(self):
		new_d, ok = QtGui.QInputDialog.getDouble(self, 'Change distance between centers', 'd=')
		
		global fish_dist_px
		fish_dist_px=new_d
		
		if ok:
			self.lbl_dist_geom.setText("distance between centers="+str(new_d))
		return new_d
	
	#Change fly ball radius
	def change_fly_radius_px(self):
		new_r, ok = QtGui.QInputDialog.getDouble(self, 'Change ball radius', 'r=')
		
		global fly_radius_px
		fly_radius_px=new_r
		
		if ok:
			self.lbl_radius_geom.setText("r="+str(new_r))
		return new_r
	
	#Change fly zinterc
	def change_fly_zinterc_px(self):
		new_z, ok = QtGui.QInputDialog.getDouble(self, 'Change z-intercept', 'z=')
		
		global fly_zinterc_px
		fly_zinterc_px=new_z
		
		if ok:
			self.lbl_radius_geom.setText("zinterc="+str(new_z))
		return new_z
	
	#Change frog radius
	def change_frog_radius_px(self):
		new_r, ok = QtGui.QInputDialog.getDouble(self, 'Change radius', 'r=')
		
		global frog_radius_px
		fly_radius_px=new_r
		
		if ok:
			self.lbl_radius_geom.setText("radius="+str(new_r))
		return new_r
	
	#Close window -> change flag
	def closeEvent(self, event):
		global geom_win_open
		geom_win_open=0
		
		
	
				
#=====================================================================================================================================
#Window for mesh parameter setting 
#=====================================================================================================================================


class win_parms_mesh(QtGui.QWidget):
	def __init__(self, parent=None):
		super(win_parms_mesh,self).__init__(parent)
	
		self.setWindowTitle('Mesh-Generator Parameter')
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		btn_change_cellvol_px=QtGui.QPushButton('Change',self)
		
		btns_mesh=[btn_change_cellvol_px]
		
		btn_change_cellvol_px.connect(btn_change_cellvol_px, QtCore.SIGNAL('clicked()'), self.change_cellvol_px)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------

		self.lbl_cellvol_px = QtGui.QLabel("Cell Volume="+str(cellvol_px), self)
		self.lbl_meshgen = QtGui.QLabel("Mesh-Generator="+mesh_gen, self)
		
		lbls_mesh=[self.lbl_cellvol_px,self.lbl_meshgen]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid_mesh = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0), (0, 1), 
			(1, 0), (1, 1)]
			
		#Defining Order of object types
		order=[1,2,1,3]
		
		#Object counters
		j=0
		k=0
		l=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
				
			if order[i]==1:
				curr_obj = lbls_mesh[k]
				k=k+1
				
			elif order[i]==2:
				curr_obj = btns_mesh[l]
				l=l+1
			
			elif order[i]==3:
				curr_obj = QtGui.QLabel('')
				
			#Adding new obj to grid layout
			grid_mesh.addWidget(curr_obj, pos[i][0], pos[i][1])
		
		self.setLayout(grid_mesh)

	
	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------
	
	
	#Change cellvol_pxume
	def change_cellvol_px(self):
		new_c, ok = QtGui.QInputDialog.getDouble(self, 'Change cell volume', 'V=',decimals=2)
		
		global cellvol_px
		cellvol_px=new_c
		
		if ok:
			self.lbl_cellvol_px.setText("Cellvolume="+str(new_c))
		return new_c
	
	def press_save_mesh(self):
		if mesh_gen=="Gmsh":
			
			#Choosing filename
			fn_save = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '.')
			
			#updating vars
			vars_mesh_names=["mesh_gen","cellvol_px"]
			vars_mesh_values=[mesh_gen,cellvol_px]
			
			#Writing into csv file
			save_sim_settings(vars_geom_names,vars_geom_values,vars_mesh_names,vars_mesh_values,vars_pde_names,vars_pde_values,vars_plot_names,vars_plot_values,fn_save)
		else:
			print "not working yet"
			
			
#=====================================================================================================================================
#Window for pde parameter setting 
#=====================================================================================================================================


class win_parms_pde(QtGui.QWidget):
	def __init__(self, parent=None):
		super(win_parms_pde,self).__init__(parent)
	
		self.setWindowTitle('PDE Parameter')
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
	
		btn_change_dt=QtGui.QPushButton('Change',self)
		btn_change_steps=QtGui.QPushButton('Change',self)
		btn_change_D=QtGui.QPushButton('Change',self)
		btn_change_k_prod=QtGui.QPushButton('Change',self)
		
		
		btns_pde=[btn_change_dt,btn_change_steps,btn_change_D,btn_change_k_prod]
		
		btn_change_dt.connect(btn_change_dt, QtCore.SIGNAL('clicked()'), self.change_dt)
		btn_change_steps.connect(btn_change_steps, QtCore.SIGNAL('clicked()'), self.change_steps)
		btn_change_D.connect(btn_change_D, QtCore.SIGNAL('clicked()'), self.change_D)
		btn_change_k_prod.connect(btn_change_k_prod, QtCore.SIGNAL('clicked()'), self.change_k_prod)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------

		self.lbl_dt = QtGui.QLabel("dt="+str(dt), self)
		self.lbl_steps = QtGui.QLabel("steps="+str(steps), self)
		self.lbl_D_nd= QtGui.QLabel("D="+str(D), self)
		self.lbl_k_prod = QtGui.QLabel("k="+str(k_prod), self)
		self.lbl_solver = QtGui.QLabel("Solver="+solver, self)
		
		
		lbls_pde=[self.lbl_dt,self.lbl_steps,self.lbl_D,self.lbl_k_prod,self.lbl_solver]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid_pde = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0), (0, 1), 
			(1, 0), (1, 1),
			(2, 0), (2, 1),
			(3, 0), (3, 1),
			(4, 0), (4, 1)]
			
		#Defining Order of object types
		order=[1,2,1,2,1,2,1,2,1,3]
		
		#Object counters
		j=0
		k=0
		l=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
			if order[i]==1:
				curr_obj = lbls_pde[k]
				k=k+1
				
			elif order[i]==2:
				curr_obj = btns_pde[l]
				l=l+1
				
			elif order[i]==3:
				curr_obj = QtGui.QLabel('')
					
			
			#Adding new obj to grid layout
			grid_pde.addWidget(curr_obj, pos[i][0], pos[i][1])
		
		self.setLayout(grid_pde)

	
	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------
	
	
	#Change dt
	def change_dt(self):
		new_dt, ok = QtGui.QInputDialog.getDouble(self, 'Change size of timesteps', 'dt=',decimals=2)
		
		global dt
		dt=new_dt
		
		if ok:
			self.lbl_dt.setText("dt="+str(new_dt))
		return new_dt

	#Change steps
	def change_steps(self):
		new_steps, ok = QtGui.QInputDialog.getInt(self, 'Change numbers of timesteps', 'steps=')
		
		global steps
		steps=new_steps
		
		if ok:
			self.lbl_steps.setText("steps="+str(new_steps))
		return new_steps

	#Change Dusivity
	def change_D(self):
		new_D, ok = QtGui.QInputDialog.getDouble(self, 'Change Dusivity', 'D=')
		
		global D
		D=new_D
		
		if ok:
			self.lbl_D.setText("D="+str(new_D))
		return new_D
	
	#Change k_prodion rate
	def change_k_prod(self):
		new_k_prod, ok = QtGui.QInputDialog.getDouble(self, 'Change k_prodion rate', 'k=')
		
		global k_prod
		k_prod=new_k_prod
		
		if ok:
			self.lbl_k_prod.setText("k="+str(new_k_prod))
		return new_k_prod
	
	#Save button
	def press_save_pde(self):
		
		#Choosing filename
		fn_save = QtGui.QFileDialog.getSaveFileName(self, 'Save File', '.')

		#Updating vars
		vars_pde_names=["solver","steps","dt","D","k_prod"]
		vars_pde_values=[solver,steps,dt,D,k_prod]
		
		#Writing into csv file
		save_sim_settings(vars_geom_names,vars_geom_values,vars_mesh_names,vars_mesh_values,vars_pde_names,vars_pde_values,vars_plot_names,vars_plot_values,fn_save)

#=====================================================================================================================================
#Window for pde parameter setting 
#=====================================================================================================================================


class win_parms_plot(QtGui.QWidget):
	def __init__(self, parent=None):
		super(win_parms_plot,self).__init__(parent)
	
		self.setWindowTitle('Plot Settings')
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		btn_add_slice=QtGui.QPushButton('Add Slice',self)
		btn_del_slice=QtGui.QPushButton('Delete Slice',self)

		btns_plot=[btn_add_slice,btn_del_slice]
		
		btn_add_slice.connect(btn_add_slice, QtCore.SIGNAL('clicked()'), self.add_slice)
		btn_del_slice.connect(btn_del_slice, QtCore.SIGNAL('clicked()'), self.del_slice)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		str_slices=''
		for i in range(slice_number):
			if i==0:
				str_slices=str_slices+str(slice_heights_px[i])
			else:
				str_slices=str_slices+','+str(slice_heights_px[i])
		

		self.lbl_nslices = QtGui.QLabel("number of slices="+str(slice_number), self)
		self.lbl_hslices = QtGui.QLabel("slices heights="+str_slices, self)
		self.lbl_output = QtGui.QLabel("Output", self)
		self.lbl_show = QtGui.QLabel("Show", self)
		self.lbl_wire = QtGui.QLabel("Wire", self)
		self.lbl_surf = QtGui.QLabel("Surface", self)
		self.lbl_conc = QtGui.QLabel("Concentrations", self)
		self.lbl_cont = QtGui.QLabel("Contour", self)
		self.lbl_all = QtGui.QLabel("All", self)
		self.lbl_compare = QtGui.QLabel("Compare with 2D", self)
		
		lbls_plot=[self.lbl_nslices,self.lbl_hslices,self.lbl_show,self.lbl_output,self.lbl_wire,self.lbl_surf,self.lbl_conc,self.lbl_cont,self.lbl_all,self.lbl_compare]
	
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------

		self.cb_show_wire = QtGui.QCheckBox('', self)
		self.cb_show_surf = QtGui.QCheckBox('', self)
		self.cb_show_conc = QtGui.QCheckBox('', self)
		self.cb_show_cont = QtGui.QCheckBox('', self)
		self.cb_show_all = QtGui.QCheckBox('', self)
		
		self.cb_out_wire = QtGui.QCheckBox('', self)
		self.cb_out_surf = QtGui.QCheckBox('', self)
		self.cb_out_conc = QtGui.QCheckBox('', self)
		self.cb_out_cont = QtGui.QCheckBox('', self)
		self.cb_out_all = QtGui.QCheckBox('', self)
		
		self.cb_compare = QtGui.QCheckBox('', self)
		
		#Fix checks on initialization
		if show_conc==1:
			self.cb_show_conc.setCheckState(QtCore.Qt.Checked)
		elif show_conc==0:
			self.cb_show_conc.setCheckState(QtCore.Qt.Unchecked)
		if show_wire==1:
			self.cb_show_wire.setCheckState(QtCore.Qt.Checked)
		elif show_wire==0:
			self.cb_show_wire.setCheckState(QtCore.Qt.Unchecked)
		if show_surf==1:
			self.cb_show_surf.setCheckState(QtCore.Qt.Checked)
		elif show_surf==0:
			self.cb_show_surf.setCheckState(QtCore.Qt.Unchecked)
		if show_cont==1:
			self.cb_show_cont.setCheckState(QtCore.Qt.Checked)
		elif show_cont==0:
			self.cb_show_cont.setCheckState(QtCore.Qt.Unchecked)
		if show_all==1:
			self.cb_show_all.setCheckState(QtCore.Qt.Checked)
		elif show_all==0:
			self.cb_show_all.setCheckState(QtCore.Qt.Unchecked)
		if out_conc==1:
			self.cb_out_conc.setCheckState(QtCore.Qt.Checked)
		elif out_conc==0:
			self.cb_out_conc.setCheckState(QtCore.Qt.Unchecked)
		if out_wire==1:
			self.cb_out_wire.setCheckState(QtCore.Qt.Checked)
		elif out_wire==0:
			self.cb_out_wire.setCheckState(QtCore.Qt.Unchecked)
		if out_surf==1:
			self.cb_out_surf.setCheckState(QtCore.Qt.Checked)
		elif out_surf==0:
			self.cb_out_surf.setCheckState(QtCore.Qt.Unchecked)
		if out_cont==1:
			self.cb_out_cont.setCheckState(QtCore.Qt.Checked)
		elif out_cont==0:
			self.cb_out_cont.setCheckState(QtCore.Qt.Unchecked)
		if out_all==1:
			self.cb_out_all.setCheckState(QtCore.Qt.Checked)
		elif out_all==0:
			self.cb_out_all.setCheckState(QtCore.Qt.Unchecked)
		if comp_flag==1:
			self.cb_compare.setCheckState(QtCore.Qt.Checked)
		elif comp_flag==0:
			self.cb_compare.setCheckState(QtCore.Qt.Unchecked)
		
		cbs_plot=[self.cb_show_wire,self.cb_out_wire,self.cb_show_surf,self.cb_out_surf,self.cb_show_conc,self.cb_out_conc,self.cb_show_cont,self.cb_out_cont,self.cb_show_all,self.cb_out_all,self.cb_compare]
		
		#Cb actions
		self.connect(self.cb_show_wire, QtCore.SIGNAL('stateChanged(int)'), self.check_show_wire)
		self.connect(self.cb_out_wire, QtCore.SIGNAL('stateChanged(int)'), self.check_out_wire)
		self.connect(self.cb_show_surf, QtCore.SIGNAL('stateChanged(int)'), self.check_show_surf)
		self.connect(self.cb_out_surf, QtCore.SIGNAL('stateChanged(int)'), self.check_out_surf)
		self.connect(self.cb_show_conc, QtCore.SIGNAL('stateChanged(int)'), self.check_show_conc)
		self.connect(self.cb_out_conc, QtCore.SIGNAL('stateChanged(int)'), self.check_out_conc)
		self.connect(self.cb_show_cont, QtCore.SIGNAL('stateChanged(int)'), self.check_show_cont)
		self.connect(self.cb_out_cont, QtCore.SIGNAL('stateChanged(int)'), self.check_out_cont)
		self.connect(self.cb_show_all, QtCore.SIGNAL('stateChanged(int)'), self.check_show_all)
		self.connect(self.cb_out_all, QtCore.SIGNAL('stateChanged(int)'), self.check_out_all)
		self.connect(self.cb_compare, QtCore.SIGNAL('stateChanged(int)'), self.check_compare)

		
		#-------------------------------------------------------------------------------------------------------------------
		#Grid
		#-------------------------------------------------------------------------------------------------------------------
		
		#Creating Grid
		grid_plot = QtGui.QGridLayout()
		
		#Defining Grid Positions
		pos = [(0, 0), (0, 1), (0, 2), 
			(1, 0), (1, 1), (1, 2),
			(2, 0), (2, 1), (2, 2),
			(3, 0), (3, 1), (3, 2),
			(4, 0), (4, 1), (4, 2),
			(5, 0), (5, 1), (5, 2),
			(6, 0), (6, 1), (6, 2),
			(7, 0), (7, 1), (7, 2),
			(8, 0), (8, 1), (8, 2)]
			
		#Defining Order of object types
		order=[1,0,2,1,0,2,0,1,1,1,3,3,1,3,3,1,3,3,1,3,3,1,3,3,1,3,0]
		
		#Object counters
		j=0
		k=0
		l=0
		m=0
		n=0
		
		for i in range(shape(pos)[0]):
			
			#Defining what next obj is
			if order[i]==1:
				curr_obj = lbls_plot[k]
				k=k+1
				
				
			elif order[i]==2:
				curr_obj = btns_plot[l]
				l=l+1
			
			elif order[i]==0:
				curr_obj= QtGui.QLabel('')
				m=m+1
			
			elif order[i]==3:
				curr_obj= cbs_plot[n]
				n=n+1
		
			#Adding new obj to grid layout
			grid_plot.addWidget(curr_obj, pos[i][0], pos[i][1])
		
		self.setLayout(grid_plot)

	
	#-------------------------------------------------------------------------------------------------------------------
	#Functions
	#-------------------------------------------------------------------------------------------------------------------
	
	#Add Slice
	def add_slice(self):
		new_slice, ok = QtGui.QInputDialog.getDouble(self, 'Enter slice height', 'h=')
		
		global slice_heights_px
		slice_heights_px.append(new_slice)
		
		global slice_number
		slice_number=slice_number+1
		
		if ok:
			str_slices=''
			for i in range(slice_number):
				if i==0:
					str_slices=str_slices+str(slice_heights_px[i])
				else:
					str_slices=str_slices+','+str(slice_heights_px[i])
			
			self.lbl_nslices.setText("n="+str(slice_number))
			self.lbl_hslices.setText("h="+str_slices)
			
		return new_slice
	
	
	#Delete Slice
	def del_slice(self):
		rem_slice, ok = QtGui.QInputDialog.getDouble(self, 'Remove slice', 'h=')
		
		global slice_heights_px
		slice_heights_px.remove(rem_slice)
		
		global slice_number
		slice_number=slice_number-1
		
		if ok:
			str_slices=''
			for i in range(slice_number):
				
				if i==0:
					str_slices=str_slices+str(slice_heights_px[i])
				else:
					str_slices=str_slices+','+str(slice_heights_px[i])
			
			self.lbl_nslices.setText("n="+str(slice_number))
			self.lbl_hslices.setText("h="+str_slices)
			
		return rem_slice
	
	#Checkbox Wire Show
	def check_show_wire(self, value):
		global show_wire
		if self.cb_show_wire.isChecked():
			show_wire=1
		else:
			show_wire=0
	
	#Checkbox Wire Out
	def check_out_wire(self, value):
		global out_wire
		if self.cb_out_wire.isChecked():
			out_wire=1
		else:
			out_wire=0
	
	#Checkbox surf Show
	def check_show_surf(self, value):
		global show_surf
		if self.cb_show_surf.isChecked():
			show_surf=1
		else:
			show_surf=0
	
	#Checkbox surf Out
	def check_out_surf(self, value):
		global out_surf
		if self.cb_out_surf.isChecked():
			out_surf=1
		else:
			out_surf=0
	#Checkbox conc Show
	def check_show_conc(self, value):
		global show_conc
		if self.cb_show_conc.isChecked():
			show_conc=1
		else:
			show_conc=0
	
	#Checkbox conc Out
	def check_out_conc(self, value):
		global out_conc
		if self.cb_out_conc.isChecked():
			out_conc=1
		else:
			out_conc=0
			
	#Checkbox cont Show
	def check_show_cont(self, value):
		global show_cont
		if self.cb_show_cont.isChecked():
			show_cont=1
		else:
			show_cont=0
	
	#Checkbox cont Out
	def check_out_cont(self, value):
		global out_cont
		if self.cb_out_cont.isChecked():
			out_cont=1
		else:
			out_cont=0
			
	#Checkbox all Show
	def check_show_all(self, value):
		global show_all
		if self.cb_show_all.isChecked():
			show_all=1
		else:
			show_all=0
	
	#Checkbox Wire Out
	def check_out_all(self, value):
		global out_all
		if self.cb_out_all.isChecked():
			out_all=1
		else:
			out_all=0
	
	#Checkbox Wire Out
	def check_compare(self, value):
		global comp_flag
		if self.cb_compare.isChecked():
			comp_flag=1
		else:
			comp_flag=0
		
def main():
		    
	#Creating application

	app = QtGui.QApplication(sys.argv)
	main_win = pyfrp_Main()
	main_win.show()
	
	
	sys.exit(app.exec_())
	
	
if __name__ == '__main__':
	main()
	

