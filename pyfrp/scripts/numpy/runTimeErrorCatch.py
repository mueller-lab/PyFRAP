import numpy as np

#np.seterr(all='raise')
with np.errstate(divide='raise'):
	try:
		np.array([1])/1.
	except FloatingPointError:
		print "bla"
