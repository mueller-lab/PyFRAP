#!/bin/bash   

#Set PYTHONUSERBASE
export PYTHONUSERBASE=/fml/ag-mueller_share/PyFRAP_and_IdealFRAPData/python-packages/
echo 'PYTHONUSERBASE:' $PYTHONUSERBASE

#Set working directory
cd /fml/ag-mueller_share/PyFRAP_and_IdealFRAPData/code
echo 'working directory:'
pwd

#Run whatever python script
python analyzeEmbryoOnCluster.py ../data/Gary/Fscin40kDa_500nM/embryoFiles/Fscin40kDa_500nM_nobeads_n0_f0_b0_g0_m0/20150818_03.emb 

#python test_final.py
#python testgmsh.py
#python hello.py
#python testRunGmshOnCluster.py
#gmsh -3 meshfiles/dome.geo
#python testMeshWithSubOnCluster.py
#python testMeshOnCluster.py
#python testCustomMeshClass.py
#python testEmbryoOnCluster.py
echo 'Done'
