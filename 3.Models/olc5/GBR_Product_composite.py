import rsgislib

import rsgislib.imagecalc

inputImages = ['20191006_Depth_Estimates.tif', '20200119_Depth_Estimates.tif', '20200213_Depth_Estimates.tif', '20201030_Depth_Estimates.tif', '20210324_Depth_Estimates.tif', 'GBR_20161115_Depth_Estimates.tif', 'GBR_20210202_Depth_Estimates.tif', 'GBR_20181125_Depth_Estimates.tif', 'GBR_20191110_Depth_Estimates.tif', 'GBR_20210513_Depth_Estimates.tif', 'GBR_20210722_Depth_Estimates.tif']

outputImage = 'GBR_comp_MAX.tif'

rsgislib.imagecalc.calcMultiImgBandStats(inputImages, outputImage, rsgislib.SUMTYPE_MAX,'GTIFF', rsgislib.TYPE_32FLOAT, 0, False)

#SUMTYPE_MEAN
#SUMTYPE_STDDEV
#SUMTYPE_MIN
