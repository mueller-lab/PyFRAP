#where_test

import numpy as np
import time

def bla(x,y):
	
	print x
	
	return x>y

#Returns true if coord is in poly
def checkInsidePoly(x,y,poly):
	
	#Taken from http://www.ariel.com.au/a/python-point-int-poly.html
	n = len(poly)
	inside =False
	
	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y
	
	return inside


tstart=time.clock()

a = 1+np.random.random(1000)*512
b = 1+np.random.random(1000)*512

print a[0:3]
print b[0:3]



c = zip(a,b)



print np.shape(c)

raw_input()

corners=[(100,100),(200,100),(200,200),(100,200),(150,150)]
idxs=[]
for i in range(len(a)):

	if checkInsidePoly(a[i],b[i],corners):
		idxs.append(i)

print len(idxs)

print time.clock()-tstart
