#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyFRAP: A Python based FRAP analysis tool box
"""

#Numpy/Scipy
#import numpy as np

#Matplotlib
#import matplotlib.pyplot as plt

#Basic PyFRAP modules
from . import modules

#PyFRAP classes
from . import subclasses

#PyFRAP GUI classes

#Only import if not RTD or not selected. RTD is currently having problems 
#with importing GUI classes. Will need to have fix for this 
#at some point

import os
import sys

if os.environ.get('READTHEDOCS', None) == 'True':
	pass
else:
	#from . import gui
	from .gui.pyfrp_app import main
	
	
	
__version__ = '1.1'
__author__ = u"Alexander Blaessle"
__license__ = "GNU GPL v3"