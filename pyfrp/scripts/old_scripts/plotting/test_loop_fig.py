import matplotlib.pyplot as plt
import numpy as np
import time

#Creating figure
fig=plt.figure()
ax=fig.add_subplot(111)
fig.show()

x=np.arange(0,30,1)

#plt.ion()
#plt.hold(True)
for i in range(20):
	
	ax.cla()
	
	y=np.random.rand(len(x))
	ax.plot(x,y,'r-')
	plt.draw()
	time.sleep(0.1)
	plt.pause(0.0001)
	

print "done"
	