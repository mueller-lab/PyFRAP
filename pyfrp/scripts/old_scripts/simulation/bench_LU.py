from scipy import linalg
import numpy as np
import time
from pysparse import superlu

tstart=time.time()
for i in range(1000000):
	A=np.random.rand(3,3)
	P,L,U=linalg.lu(A)
	
	
	
print "Took", time.time()-tstart

