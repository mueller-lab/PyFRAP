# PyFRAP
PyFRAP: A Python based FRAP analysis tool box

## Features

- Import FRAP datasets from timelapse experiments and analyze image data with various options such as
	+ various filters
	+ background substraction
	+ illumination correction
- Simulate the FRAP experiment with exact interpolated initial conditions
- Fit simulated experiment to analyzed data and extract diffusion coefficient
- Statistical analysis of fitting results
- Hierarchical data structure making data exchange/sharing easy

# Getting Started

## Installation

If you have Python2.7 already installed, simply 

	git clone https://github.com/alexblaessle/PyFRAP
	
and install via:

	python setup.py install
	
For a full installation documentation, have a look at the documentation.

### Requirements

PyFRAP depends on 

- numpy>=1.8.2
- matplotlib>=1.4.3
- scipy>=0.13.3
- scikit-image>=0.11.3
- fipy>=3.1
- PyQT4>4.10.4
- vtk>=5.8.0
- colorama>=0.2.5
- meshio>=1.2.1
- pickle
- gmsh (compiled with TetGen Algorithm)

and a bunch of standard Python2.7 libraries such as

- time
- platform
- sys
- shutil
- tempfile
- gc
- csv

## API

The API of PyFRAP can be found [here](http://pyfrp.readthedocs.org/en/latest/pyrw.html#submodules "toAPI") .

## Documentation

<!-- The Documentation of pyrw can be found [here](http://pyrw.readthedocs.org/en/latest/pyrw.html "toAPI") . -->

## Example

