import numpy as np
import cv2
from cv2 import cv
import scipy
from scipy import ndimage

fn='/home/alex_loc/Documents/Research/PyFRAP/Ideal_FRAP/zStackData/20160119_zstack_nobeads/01/No_beads_010017.tif'

im_col = cv2.imread(fn)
im = cv2.imread(fn,cv2.CV_LOAD_IMAGE_GRAYSCALE)

im = np.where(im>100,0,255).astype(np.uint8)
im = cv2.erode(im, None,iterations=8)
im_label, num = ndimage.label(im)



for label in xrange(1, num+1):
	points = np.array(np.where(im_label==label)[::-1]).T.reshape(-1,1,2).copy()
	
	#print points.shape
	#print points[:10,:,:]
	#raw_input()
	
	for i in range(points.shape[0]):
		center=tuple(list(points[i,:,:][0]))
		print center
		cv2.circle(im_col,center,1,(0,255,0),-1)
		
	center,radius= cv2.minEnclosingCircle(points)
	center=tuple(np.asarray(center).astype(np.int))
	radius=int(radius)
	#lines = np.array(cv2.cv.BoxPoints(rect)).astype(np.int)
	#if any([np.any(lines[:,0]<=0), np.any(lines[:,0]>=im.shape[1]-1), np.any(lines[:,1]<=0), np.any(lines[:,1]>=im.shape[0]-1)]):
		#cv2.drawContours(im_col,[lines],0,(0,0,255),1)
	#else:
		#cv2.drawContours(im_col,[lines],0,(255,0,0),1)
	print im_col
	print center,radius
	#cv2.circle(im_col,center,radius,(255,0,0),-1)

print im_col

cv2.imshow('im',im_col)
#cv2.imwrite('rects.png',im_col)
cv2.waitKey()