import numpy as np

x=np.random.random(1000)*10
y=np.random.random(1000)*10

offset=[2,2]
sidelength=5

upper=offset[0]+sidelength

indSquare=np.where((offset[0]<x)  & (x<upper) & (offset[1]<y) & (y<offset[1]+sidelength))[0]

print len(indSquare)