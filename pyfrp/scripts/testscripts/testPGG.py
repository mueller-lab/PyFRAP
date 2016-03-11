###TestScript for pyfrp_gmsh_geometry

import pyfrp_gmsh_geometry as pgg 

d=pgg.domain()


#Define Vertices of upper circle
centerUp=d.addVertex([0.,0.,0.])

rightUp=d.addVertex([1.,0.,0.])
topUp=d.addVertex([0.,1.,0.])
leftUp=d.addVertex([-1.,0.,0.])
bottomUp=d.addVertex([0.,-1.,0.])

#Define Vertices of lower circle
centerLow=d.addVertex([0.,0.,-2.])

rightLow=d.addVertex([1.,0.,-2.])
topLow=d.addVertex([0.,1.,-2.])
leftLow=d.addVertex([-1.,0.,-2.])
bottomLow=d.addVertex([0.,-1.,-2.])

#Define Lines connecting circles
d.addLine(rightUp,rightLow)
d.addLine(topUp,topLow)
d.addLine(leftUp,leftLow)
d.addLine(bottomUp,bottomLow)

#Define Arcs of upper circle
d.addArc(rightUp,centerUp,topUp)
d.addArc(topUp,centerUp,leftUp)
d.addArc(leftUp,centerUp,bottomUp)
d.addArc(bottomUp,centerUp,rightUp)

#Define Arcs of lower circle
d.addArc(rightLow,centerLow,topLow)
d.addArc(topLow,centerLow,leftLow)
d.addArc(leftLow,centerLow,bottomLow)
d.addArc(bottomLow,centerLow,rightLow)

d.draw()



raw_input()

