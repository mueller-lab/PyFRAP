import sys
import cv2
import numpy
import random
from scipy.ndimage import label

def segment_on_dt(img):
	dt = cv2.distanceTransform(img, 2, 3) # L2 norm, 3x3 mask
	dt = ((dt - dt.min()) / (dt.max() - dt.min()) * 255).astype(numpy.uint8)
	dt = cv2.threshold(dt, 100, 255, cv2.THRESH_BINARY)[1]
	lbl, ncc = label(dt)

	lbl[img == 0] = lbl.max() + 1
	lbl = lbl.astype(numpy.int32)
	cv2.watershed(cv2.cvtColor(img, cv2.COLOR_GRAY2BGR), lbl)
	lbl[lbl == -1] = 0
	return lbl


img = cv2.cvtColor(cv2.imread(sys.argv[1]), cv2.COLOR_BGR2GRAY)
img = cv2.threshold(img, 0, 255, cv2.THRESH_OTSU)[1]
img = 255 - img # White: objects; Black: background

ws_result = segment_on_dt(img)
# Colorize
height, width = ws_result.shape
ws_color = numpy.zeros((height, width, 3), dtype=numpy.uint8)
lbl, ncc = label(ws_result)
for l in xrange(1, ncc + 1):
	a, b = numpy.nonzero(lbl == l)
	if img[a[0], b[0]] == 0: # Do not color background.
		continue
	rgb = [random.randint(0, 255) for _ in xrange(3)]
	ws_color[lbl == l] = tuple(rgb)

cv2.imwrite(sys.argv[2], ws_color)