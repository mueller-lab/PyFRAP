#! /bin/bash

# This script installs all dependencies of PyFRAP for OSX, assuming Anaconda is installed.

# Make folder
# mkdir $1
# cd $1

# Downgrade to PyQt4
conda remove -y qt pyqt
conda install -y qt=4 pyqt=4.10.4

# Reinstall MPL
conda remove -y pywavelets
conda install -y matplotlib=1.4.3

# Reinstall skimage
conda install -y scikit-image=0.11.3

# Install extra packages
conda install -y vtk
conda install -y colorama
conda install -y -c guyer fipy

# Install extra packages using 
pip install numpy-stl
pip install solidpython
pip install wget

# 
#git clone https://github.com/alexblaessle/PyFRAP.git
#cd PyFRAP-master

python setup.py install --gmsh --fiji