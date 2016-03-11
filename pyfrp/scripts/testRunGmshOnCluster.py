###Script to debug runGmsh function on cluster

import pyfrp_gmsh_module

fnGeo='meshfiles/cylinder.geo'
fnMsh='meshfiles/cylinder.msh'
debug=True


pyfrp_gmsh_module.runGmsh(fn,fnOut=fnMsh,debug=debug)