#!/bin/bash   

#Define folder
folder=$1

echo "Analyzing all embryo files in folder" $folder
	
for d in $folder/*;
do
	
	echo "=================================================="
	echo "Processing subfolder" $d
	
	for f in $d/*.emb;
	do
		
		echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		
		#Convert filepath to absolute path
		fileSub=$(readlink -f $f)
		
		echo "file = "  $fileSub
		
		#Generate folders for sterr and stoud
		mkdir -p $d/sterr
		mkdir -p $d/stout
# 		
		#Grab filename base
		outname=$(basename $f)
		
		#Define paths for sterr/stout
		sterrPath=$(readlink -f "$d/sterr/$outname.sterr")
		stoutPath=$(readlink -f "$d/stout/$outname.stout")
		
	
		#Delete old stout/sterr files if they exist
		rm $sterrPath
		rm $stoutPath
		
		#Print out stoutpath
		echo "stout = $stoutPath"
		echo "sterr = $sterrPath"
		 
		#Submit
		qsub -l h_vmem=8G -l h_rt=00:60:00 -o $stoutPath -e $sterrPath analyzeEmbryoOnCluster.sh $fileSub
		
	done
	
# 	
# 	
done

