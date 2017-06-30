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
- Comprehensive GUI, making almost all PyFRAP tools available

## Installation

If you have Python2.7 already installed, simply 

	git clone https://github.com/alexblaessle/PyFRAP
	
and install via:

	python setup.py install --user
	
We highly recommend installing with the *--user* option, since PyFRAP needs to read/write data files in the installation folder. In some cases, this
might lead to file permission issues.
For a full installation documentation, have a look at the [wiki](https://github.com/alexblaessle/PyFRAP/wiki/Installation).

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
- wget>=3.2
- gmsh (compiled with TetGen Algorithm) MUST BE Version 2.14.0!

## Getting Started

### Running PyFRAP GUI

PyFRAP comes with comprehensive GUI. To start the GUI, simply go to pyfrp/ and doubleclick on 

	PyFRAP.py
	
or run

	python pyfrp/PyFRAP.py

If you are already in a python session, you can simply run 

	import pyfrp
	pyfrp.main()

Note that in the latter method PyFRAP's stdout might get redirected to the python shell you are executing it from.

### Running PyFRAP from the command line

```python
	import pyfrp
```

### Using PyFRAP GUI to analyze a FRAP experiments

Check out the PyFRAP wiki's [First Steps Section](https://github.com/alexblaessle/PyFRAP/wiki/FirstSteps).

## API

The API of PyFRAP can be found [here](http://pyfrap.readthedocs.org/en/latest/) .

## Documentation

To learn more about PyFRAP, check out the PyFRAP [wiki](https://github.com/alexblaessle/PyFRAP/wiki)





