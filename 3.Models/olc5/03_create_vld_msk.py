import rsgislib.imageutils
import os

os.environ["RSGISLIB_IMG_CRT_OPTS_GTIFF"] = "TILED=YES:COMPRESS=LZW"

input_img = "GBR_20201030.kea"
out_vld_img = "Study_site_vmsk.tif"

rsgislib.imageutils.genValidMask(input_img, out_vld_img, 'GTIFF', 0)

