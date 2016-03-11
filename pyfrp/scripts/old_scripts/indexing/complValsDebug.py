import numpy as np
import pyfrp_misc_module
import time

#Generate random samples
s=40000

a=np.random.randint(10000,size=s)
b=np.random.randint(10000,size=s)

tStart=time.clock()

a1=pyfrp_misc_module.complValsSimple(a,b)

print "Simple done after: ", time.clock()-tStart, " removing ", len(a)-len(a1), " matches."

tStart=time.clock()

a2=pyfrp_misc_module.complValsFast(a,b)

print "Fast done after: ", time.clock()-tStart, " removing ", len(a)-len(a2), " matches."