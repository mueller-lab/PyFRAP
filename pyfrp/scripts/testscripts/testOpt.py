from scipy.optimize import fmin
from scipy.optimize import minimize

def func(x,a,b):
	return (x[0]-b)**2+(x[1]-a)**2



print fmin(func,[1,1],args=(2,3,))

print minimize(func,[1,1],args=(2,3,))
