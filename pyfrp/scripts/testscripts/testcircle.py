import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def circlefit3d(p1,p2,p3):

	# Default values
	center = []
	rad = 0
	v1n=[]
	v2nb=[]

	# Start calculation
	# v1, v2 describe the vectors from p1 to p2 and p3, resp.

	v1 = p2 - p1
	v2 = p3 - p1

	# l1, l2 describe the lengths of those vectors
	l1 = np.sqrt((v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]))
	l2 = np.sqrt((v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]))
	
	if l1==0 or l2==0: ##ok<OR2>
		print 'Error using ==> cirlefit3d. Corresponding input points must not be identical.'
		rad = -4
		return center,rad,v1n,v2nb
	
	# v1n, v2n describe the normalized vectors v1 and v2
	v1n = v1/l1
	v2n = v2/l2
	
	# nv describes the normal vector on the plane of the circle
	nv = np.array([v1n[1]*v2n[2] - v1n[2]*v2n[1] , v1n[2]*v2n[0] - v1n[0]*v2n[2] , v1n[0]*v2n[1] - v1n[1]*v2n[0]])
	
	if sum(abs(nv))<1e-5:
		print "Warning using ==> cirlefit3d. Some corresponding input points are nearly collinear."
	
	# v2nb: orthogonalization of v2n against v1n
	dotp = v2n[0]*v1n[0] + v2n[1]*v1n[1] + v2n[2]*v1n[2]
	v2nb = v2n-dotp*v1n
	
	# normalize v2nb
	l2nb = np.sqrt((v2nb[0]*v2nb[0]+v2nb[1]*v2nb[1]+v2nb[2]*v2nb[2]))
	v2nb = v2nb/l2nb

	# remark: the circle plane will now be discretized as follows
	#
	# origin: p1                    normal vector on plane: nv
	# first coordinate vector: v1n  second coordinate vector: v2nb

	# calculate 2d coordinates of points in each plane
	# p1_2d = zeros(n,2); # set per construction
	# p2_2d = zeros(n,2);p2_2d[0] = l1; # set per construction
	p3_2d = np.zeros(np.shape(p1)) # has to be calculated
	for i in range(3):
		p3_2d[0] = p3_2d[0] + v2[i]*v1n[i]
		p3_2d[1] = p3_2d[1] + v2[i]*v2nb[i]
	
	# calculate the fitting circle 
	# due to the special construction of the 2d system this boils down to solving
	# q1 = [0,0], q2 = [a,0], q3 = [b,c] (points on 2d circle)
	# crossing perpendicular bisectors, s and t running indices:
	# solve [a/2,s] = [b/2 + c*t, c/2 - b*t]
	# solution t = (a-b)/(2*c)

	a = l1
	b = p3_2d[0]
	c = p3_2d[1]
	t = 0.5*(a-b)/c
	scale1 = b/2 + c*t
	scale2 = c/2 - b*t
	
	# centers
	center = p1 + scale1*v1n + scale2*v2nb
	
	# radii
	rad = np.sqrt((center[0]-p1[0])**2+(center[1]-p1[1])**2+(center[2]-p1[2])**2)

	return center,rad,v1n,v2nb

def getNormVector(vec1,vec2):
	return np.cross(vec1,vec2)


def intersectPlane(p,n,r0,nPlane):
	
	if np.dot(n,nPlane)==0 or np.dot(n,nPlane)==1:
		return None
	
	lambd=(np.dot(nPlane*r0)-np.dot(nPlane,p))/(np.dot(n,nPlane))
	
	return p+lambd*n

def getOrthProj(p,n):
	
	orig=np.array([0.,0.,0.])
	vec1=np.array([1.,0.,0.])
	vec2=np.array([0.,1.,0.])
	
	nPlane=getNormVector(vec1,vec2)
	
	inter=intersectPlane(p,n,orig,nPlane)
	
	return inter

def turnCircle(center,start,end):
	
	radius=np.linalg.norm(center-start)
	
	v1=start-center
	v2=end-center
	
	v1n = v1/np.linalg.norm(v1)
	v2n = v2/np.linalg.norm(v2)
	
	dotp = v2n[0]*v1n[0] + v2n[1]*v1n[1] + v2n[2]*v1n[2]
	v2nb = v2n-np.dot(v2n,v1n)*v1n
	
	pOffset=radius*v2nb
	
	angleOffset=getAngle(pOffset,v1)
	angle=getAngle(v1,v2)
	
	a = np.linspace(angleOffset,angleOffset+angle,1000)
	
	x = center[0]+np.sin(a)*radius*v1n[0]+np.cos(a)*radius*v2nb[0]
	y = center[1]+np.sin(a)*radius*v1n[1]+np.cos(a)*radius*v2nb[1]
	z = center[2]+np.sin(a)*radius*v1n[2]+np.cos(a)*radius*v2nb[2]
	
	return x,y,z
	
def getAngle(vec1,vec2):
	a=np.arccos(np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
	if a<0:
		return getAngle(vec2,vec1)
	return a
	
center=np.array([0,0,0])
start=np.array([0,1,0])
end=np.array([0,0,1])

x,y,z=turnCircle(center,start,end)

print np.shape(x)

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.plot(x,y,zs=z)

plt.show()

