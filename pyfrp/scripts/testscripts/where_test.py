import numpy as np


def checkInsideCircle2(xyvec,center,radius):
	sol=[]
	for xy in xyvec:
		if np.sqrt((xy[0]-center[0])**2+(xy[1]-center[1])**2)<radius:
			sol.append(True) 
		else:
			sol.append(False)
	return sol

def checkInsideCircle(x,y,center,radius):
	return np.sqrt((x-center[0])**2+(y-center[1])**2)<radius
	
center=[256,256]
radius=300

x=np.arange(-500,1500,1)
y=np.arange(-500,1500,1)

#xy=zip(x,y)
xy=[list(a) for a in zip(x,y)]

xy=np.asarray(xy)



#xy=np.vstack((x,y))

#print np.shape(xy)
#raw_input()

#xy=np.sqrt(np.sum(((xy.T-np.array(center))**2),axis=1))

#print np.shape(xy)

print len(np.where(checkInsideCircle2(xy,center,radius))[0])

print sum(checkInsideCircle(x,y,center,radius))

#print np.where(xy<radius)
