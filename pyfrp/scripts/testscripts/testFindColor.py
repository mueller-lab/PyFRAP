###Script for testing findColor for pyfrapterminal

def initColorKeys():
	colorKeys={}
	
	colorKeys['[32m']='green'
	colorKeys['[39m']='black'
	colorKeys['[33m']='yellow'
	colorKeys['[31m']='red'
	return colorKeys

def findColors(line,colorDic):
	
	found=0
	
	colors=[]
	strings=[]
	
	while found!=-1:
		
		fs=[]
		cs=[]
		for key in colorDic.keys(): 
			ind=line.find(key,found)
			if ind!=-1:
				fs.append(ind)
			else: 
				fs.append(str(ind))
		
		foundOld=found
		found=min(fs)
			
		if len(colors)==0 and found!='-1':
			colors.append(colorDic.values()[fs.index(found)])
		elif len(colors)>0:
			offset=len(colorDic.keys()[fs.index(found)])-1
			if found=='-1':
				strings.append(line[foundOld+offset:])
			else:
				strings.append(line[foundOld+offset:found])
				colors.append(colorDic.values()[fs.index(found)])
			
		found=int(found)
		if found!=-1:
			found=found+1
		
	return colors,strings

line=r"[31mERROR: [39m lala"

print line

colorDic=initColorKeys()

colors,strings=findColors(line,colorDic)

print colors
print strings