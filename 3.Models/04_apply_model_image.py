import joblib
import regress_sklearn
import os

os.environ["RSGISLIB_IMG_CRT_OPTS_GTIFF"] = "TILED=YES:COMPRESS=LZW"


# Model file path
regrs_mdl_file = 'IS_Ber_ET_Model.joblib'

# Load the model
regrs_mdl = joblib.load(regrs_mdl_file)

# Input Image
metrics_img = 'GBR_20201030.kea'
vld_img = 'Study_site_vmsk.tif'

# Image bands in the metrics image used by the model.
metrics_band_idxs = [1,2,3,4]

# Output Image
out_img = '20201030_Depth_Estimates.tif'
out_band_names = ['Depth']

# Apply the model to the image data
regress_sklearn.apply_regress_sklearn_mdl(regrs_mdl, 1, metrics_img, metrics_band_idxs, vld_img, 1, out_img, gdalformat='GTIFF', out_band_names=out_band_names, calc_stats=True, out_no_date_val=0.0)



