#PyFRAP Modules
import pyfrp_sim_module as pyfrp_sim
import pyfrp_misc_module as pyfrp_misc
import pyfrp_img_module as pyfrp_img
import pyfrp_fit_module as pyfrp_fit
import pyfrp_gmsh_module as pyfrp_gmsh
import pyfrp_plot_module as pyfrp_plot

from embryo import *
from molecule import *

mol=molecule("blub")
fn_load="/home/alex_loc/Documents/Research/florencia/fas-gfp.pk"
mol=mol.load_molecule(fn_load)

emb=mol.embryos[1]
emb.debug_analysis=1
emb = pyfrp_img.analyze_dataset(emb)

mol.save_molecule(fn_load)