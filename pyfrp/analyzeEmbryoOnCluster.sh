#!/bin/bash   

#Set proper PYTHONUSERBASE
export PYTHONUSERBASE=/fml/ag-mueller_share/PyFRAP_and_IdealFRAPData/python-packages/
echo 'Set PYTHONUSERBASE to:' $PYTHONUSERBASE

#Set WD
cd /fml/ag-mueller_share/PyFRAP_and_IdealFRAPData/code
echo 'Working Directory:' $(pwd)


#Print out date
echo '==============================================================='
echo 'Started job at ' $(date)
echo 'Will analyze file ' $1
echo '==============================================================='

#Execute python script
python analyzeEmbryoOnCluster.py $1
