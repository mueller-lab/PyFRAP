###Script to compare flattenDatasets

import os
import pyfrp_img_module
import pyfrp_misc_module
import pyfrp_plot_module

import matplotlib.pyplot as plt

#Some parameters
center=[256,256]

#Bookkeeping
imgs=[]
imgAxes=[]

#Define datasets
fnFlatten=['../Ideal_FRAP/flatteningData/20152011_flat/01/','../Ideal_FRAP/flatteningData/20152011_flat/02/','../Ideal_FRAP/flatteningData/20160105_flat/01/'
	   ,'../Ideal_FRAP/flatteningData/20160201_flat/01/','../Ideal_FRAP/flatteningData/20160201_flat/02/','../Ideal_FRAP/flatteningData/20160201_flat/03/','../Ideal_FRAP/backgroundData/20160105_bkgd/01/']
#fnFlatten=['../Ideal_FRAP/backgroundData/20160105_bkgd/01/']

#Make figure for images and hist
figImg,axesImg=pyfrp_plot_module.makeSubplot([len(fnFlatten),2])

#Make figure for radial distribution
figRad,axesRad=pyfrp_plot_module.makeSubplot([1,3])

#Loop through all datsets
for i,fn in enumerate(fnFlatten):
	
	#Load files
	fileList=pyfrp_misc_module.getSortedFileList(fn,'.tif')
	
	fileList=fileList[:500]
	
	#Compute mean img
	meanImg=pyfrp_img_module.computeMeanImg(fn,fileList,'uint16')
	
	#Compute flattening mask
	flatMask=pyfrp_img_module.computeFlatMask(meanImg)
	
	#Show results
	pyfrp_img_module.showImgAndHist(meanImg,axes=axesImg[2*i:2*i+2],fixSize=False,binSize=1.,binsFit=False,title=[fn,''])
	
	#Append imgs and axes for scale normalization
	imgs.append(meanImg)
	imgAxes.append(axesImg[2*i])
	
	#Plot radial distribution
	pyfrp_img_module.plotRadialProfile(flatMask,center,ax=axesRad[0],color=[i/len(fnFlatten),1./(i+1.1),(i+1)/len(fnFlatten)],linestyle='-')
	
	#Plot radial histogram
	name=fn.replace("../Ideal_FRAP/flatteningData/","")
	ax,ax2,binsMid,histY,binY = pyfrp_img_module.plotRadialHist(flatMask,center,axes=axesRad[1],color=[i/len(fnFlatten),1./(i+1.1),(i+1)/len(fnFlatten)],linestyle='-',label=name,legend=True,fullOutput=True)
	
	#Plot timeseries
	meanIntensities = pyfrp_img_module.getMeanIntensitiesImgs(fn,fileList,'uint16')
	axesRad[2].plot(meanIntensities,color=[i/len(fnFlatten),1./(i+1.1),(i+1)/len(fnFlatten)])
	plt.draw()
	
	print fn, binY[0]/binY[-1]
	#print histY[0]/histY[-1]
	
	#print np.mean(meanIntensities)
	
	
#Adjust img range for comparison
minImg,maxImg=pyfrp_img_module.getCommonRange(imgs)
pyfrp_plot_module.adjustImshowRange(imgAxes,vmin=minImg,vmax=maxImg)

#figDev,axesDev=pyfrp_plot_module.makeSubplot([1,1])
#axesDev[0].imshow(imgs[0]/imgs[1])
#plt.draw()



#Redraw everything
for ax in axesImg:
	pyfrp_plot_module.redraw(ax)

for ax in axesRad:
	pyfrp_plot_module.redraw(ax)

raw_input()	
	
	
	






