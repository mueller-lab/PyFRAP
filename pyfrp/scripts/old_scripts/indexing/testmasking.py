import numpy as np
import pyfrp_idx_module


#Define parameters
sidelength=100
offset=[100,100]
res=512

#Generate mask
idxX,idxY=pyfrp_idx_module.getSquareIdxImg(offset,sidelength,res,debug=False)

#To mask
mask=np.zeros((res,res))
mask=pyfrp_idx_module.ind2mask(mask,idxX,idxY,1)

#To bool
mask=mask.astype(bool)

#idx grid
x=np.arange(res)
y=np.arange(res)	
X,Y=np.meshgrid(x,y)

#Slice idx grid
idxX_new=X[mask].flatten()
idxY_new=Y[mask].flatten()

#Put in new mask
mask_new=np.zeros((res,res))
mask_new=pyfrp_idx_module.ind2mask(mask,idxX_new,idxY_new,1)

#Check if sum is 0
print (mask-mask_new).sum()

#Check if same length
print len(idxX), len(idxX_new)
print len(idxY), len(idxY_new)



