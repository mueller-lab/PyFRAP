### Script to do statistics of how big PyFRAP is

import os

cwd=os.getcwd()

files=os.listdir(cwd)

nModules=0
nClasses=0
nFunctions=0
nLines=0
nMethods=0
nGUI=0
nEmpty=0
nFiles=0

cwd=cwd+'/'

files.sort()

for fn in files:
	
	nLinesLast=nLines
	nEmptyLast=nEmpty
	
	
	if fn.endswith('py') and fn.startswith('pyfrp'):
		
		nFiles=nFiles+1
		
		if "gui" in fn:
			nGUI=nGUI+1
		
		if "module" in fn:
			nModules=nModules+1
			
		with open(cwd+fn,'r') as f:
			inClass=False
			
			for line in f:
				
				if len(line.strip())==0:
					nEmpty=nEmpty+1
					continue
				
				nLines=nLines+1
				
				if line.strip().startswith('class'):
					inClass=True
					nClasses=nClasses+1
				if line.strip().startswith('def'):
					if inClass:
						nMethods=nMethods+1
					else:
						nFunctions=nFunctions+1
		
		print fn, nLines-nLinesLast, nEmpty-nEmptyLast,  nLines-nLinesLast + nEmpty-nEmptyLast
		
		
print "nFiles = ",nFiles
print "nModules = ",nModules
print "nGUI = ",nGUI
print "nClasses = ",nClasses
print "nMethods = ", nMethods
print "nFunctions = ", nFunctions
print "nLines = ", nLines
print "nEmpty = ", nEmpty





				
			
		
		
		
		