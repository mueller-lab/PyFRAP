import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np

# create an ellipse
el = matplotlib.patches.Ellipse((50,-23), 10, 13.7, 30, facecolor=(1,0,0,.2), edgecolor='none')

# calculate the x and y points possibly within the ellipse
y_int = np.arange(-30, -15)
x_int = np.arange(40, 60)

# create a list of possible coordinates
g = np.meshgrid(x_int, y_int)
coords = list(zip(*(c.flat for c in g)))

# create the list of valid coordinates (from untransformed)
ellipsepoints = np.vstack([p for p in coords if el.contains_point(p, radius=0)])

# just to see if this works
fig = plt.figure()
fig.show()
ax = fig.add_subplot(111)
ax.add_artist(el)
ep = np.array(ellipsepoints)
ax.plot(ellipsepoints[:,0], ellipsepoints[:,1], 'ko')

plt.draw()
raw_input()