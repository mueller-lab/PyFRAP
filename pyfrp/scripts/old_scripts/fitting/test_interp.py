from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0, 10,1)
y = np.exp(-x/3.0)
f = interpolate.interp1d(x, y,bounds_error=False)
xnew = np.arange(0,9, 0.1)
ynew = f(xnew)   
 #use interpolation function returned by `interp1d`
plt.plot(x, y, 'o', xnew, ynew, '-')
#plt.plot(xnew, ynew, '-')

plt.show()