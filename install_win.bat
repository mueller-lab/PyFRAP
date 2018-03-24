REM This script installs all dependencies of PyFRAP for Windows, assuming Anaconda and the JDK are installed.

REM Downgrade to PyQt4
conda remove -y qt pyqt
conda install -y qt=4 pyqt=4.10.4

REM Reinstall MPL
conda remove -y pywavelets
conda install -y matplotlib=1.4.3

REM Install extra packages using 
python -m pip install --upgrade pip
pip install colorama
pip install numpy-stl
pip install solidpython
pip install wget
pip install fipy
pip install python-bioformats

REM Install PyFRAP
python setup.py install 

REM Reinstall skimage and install vtk
conda install -y scikit-image=0.11.3 vtk