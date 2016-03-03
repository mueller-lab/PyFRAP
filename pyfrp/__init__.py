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
from . import pyfrp_IO_module
from . import pyfrp_term_module
from . import pyfrp_misc_module
from . import pyfrp_plot_module
from . import pyfrp_img_module
from . import pyfrp_optimization_module
from . import pyfrp_stats_module
from . import pyfrp_fit_module
from . import pyfrp_stats_module
from . import pyfrp_gmsh_module
from . import pyfrp_gmsh_IO_module
from . import pyfrp_sim_module
from . import pyfrp_integration_module
from . import pyfrp_idx_module

#PyFRAP classes
from . import pyfrp_analysis
from . import pyfrp_conf
from . import pyfrp_fit
from . import pyfrp_mesh
from . import pyfrp_simulation
from . import pyfrp_geometry
from . import pyfrp_gmsh_geometry
from . import pyfrp_ROI
from . import pyfrp_embryo
from . import pyfrp_molecule

#Obsolete/Not-integrated modules  
#from . import pyfrp_zstack_module

__version__ = '1.0'
__author__ = u"Alexander Blaessle"
__license__ = "GNU GPL v3"