#Script testing howto grab pixels within a shape from matplotlib

import matplotlib.pyplot as plt
import matplotlib.patches as ptc

from numpy import *

import time

#Create figure
#fig=plt.figure()
#fig.show()

ts=time.clock()

#Create image
img=zeros((512,512)).T

#Define corners of polygon 
corners=[(100,100),(200,100),(200,200),(100,200),(150,150)]
corners=asarray(corners)



poly = ptc.Polygon(corners,edgecolor='r',facecolor=(1,0,0,.2),)


x_int=arange(1,img.shape[0]+1,1)
y_int=arange(1,img.shape[1]+1,1)

g = meshgrid(x_int, y_int)



coords = list(zip(*(c.flat for c in g)))

pts = vstack([p for p in coords if poly.contains_point(p, radius=0)])

print time.clock()-ts
ts=time.clock()

#ax.plot(pts[:,0], pts[:,1], 'g*')
#ax.add_patch(poly)
#ax2.plot(pts[:,0], pts[:,1], 'g*')
#ax.plot(pts[:,0], pts[:,1], 'g*')
#plt.draw()


#Create circle
circ = ptc.Circle([0.,15.],50,edgecolor='y',facecolor=(1,0,0,.2),)


pts = vstack([p for p in coords if circ.contains_point(p, radius=0)])
print time.clock()-ts
#ax.add_patch(circ)
#ax2.plot(pts[:,0], pts[:,1], 'm*')
#ax.plot(pts[:,0], pts[:,1], 'm*')
#plt.draw()


raw_input()



