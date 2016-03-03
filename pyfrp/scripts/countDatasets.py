### Script to calculate number of datasets for ideal frap

import os
import sys


def filterFolder(p,l):
	
	if p[-1]!="/":
		p=p+"/"
	
	lnew=[]
	for i in l:
		if os.path.isdir(p+i):
			lnew.append(i)	
	return lnew	

folder=sys.argv[1]

molFolders=filterFolder(folder,os.listdir(folder))


nEmbryo=0
nBeads=0
nNoBeads=0
nTotal=0


for molFolder in molFolders:
	
	
	nTotalBefore=nTotal
	
	tempFolders=filterFolder(folder+molFolder,os.listdir(folder+molFolder))
	
	for tempFolder in tempFolders:
		
		if tempFolder=='beads': 
			embFolders=filterFolder(folder+molFolder+"/"+tempFolder,os.listdir(folder+molFolder+"/"+tempFolder))
			nBeads=nBeads+len(embFolders)
		elif tempFolder=='nobeads': 
			embFolders=filterFolder(folder+molFolder+"/"+tempFolder,os.listdir(folder+molFolder+"/"+tempFolder))
			nNoBeads=nNoBeads+len(embFolders)
		elif tempFolder=='embryo': 
			embFolders=filterFolder(folder+molFolder+"/"+tempFolder,os.listdir(folder+molFolder+"/"+tempFolder))
			nEmbryo=nEmbryo+len(embFolders)
	
	nTotal=nBeads + nEmbryo + nNoBeads
	
	print "Counted ", nTotal - nTotalBefore, " in ", molFolder  
	
print "nEmbryo = ", nEmbryo
print "nNoBeads = ", nNoBeads
print "nBeads = ", nBeads

print "Total = ", nBeads + nEmbryo + nNoBeads


			




				
			
		
		
		
		