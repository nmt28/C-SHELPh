#!/usr/bin/env python
# coding: utf-8


#####
# This group of functions processes ICESAT2 data and creates a bathymetric model. 
# To do this, it follows a number of steps in the form of functions, including:
# 1. Reading data (ReadATL03())
# 2. Orthometrically correcting the dataset (OrthometricCorrection())
# 3. Pulling down the data segment ID (getAtl03SegID())
# 4. Bin the data along latitudinal and height gradients (bin_data())
# 5. Calculate sea height (get_sea_height())
# 6. Get water temperature (get_water_temp())
# 7. Correct bathymetry surface for refraction (RefractionCorrection())
# 8. Calculate bathymetry height (get_bath_height())
# 9. Produce figures (produce_figures())
#####

import numpy as np
import h5py as h5
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import matplotlib.pyplot as plt
import pyproj
from pyproj import Proj
from pyproj import Transformer
import pandas as pd
import argparse
import os
import subprocess
import time
import utm
import math
import fiona
import geopandas
import netCDF4
from datetime import datetime
# Import functions from utils.py
from utils import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Specify the input ICESAT H5 file")
    parser.add_argument("-l", "--laser", type=str, help="Specify the ICESAT-2 laser number (1, 2 or 3)")
    parser.add_argument("-th", "--thresh", type=int, help="Specify the threshold percentage")
    parser.add_argument("-o", "--output", type=str, required = False, help="Specify the output location")
    parser.add_argument("-lr", "--lat_res", type=float, default = 10, help="Specify the latitudinal resoltuion (normally 10)")
    parser.add_argument("-hr", "--h_res", type=float, default = 0.5, help="Specify the height resolution (normally 0.5)")
    parser.add_argument("-wt", "--waterTemp", type=float, default = 20, help="Specify the water temperature in degrees C")
    parser.add_argument("-slat", "--start_lat", type=float, required = False, help="Specify the start latitude")
    parser.add_argument("-elat", "--end_lat", type=float, required = False, help="Specify the stop latitude")
    
    args = parser.parse_args()
    
    if args.input == None:
        print('MISSING H5')
        os._exit(1)
    elif args.laser == None:
        print('MISSING LASER NUMBER')
        os._exit(1)
    elif args.waterTemp == None:
        print('MISSING WATER TEMP')
        os._exit(1)
    elif args.lat_res == None:
        print('MISSING LATITUDINAL RESOLUTION')
        os._exit(1)
    elif args.h_res == None:
        print('MISSING HEIGHT RESOLUTION')
        os._exit(1)
    elif args.thresh == None:
        print('MISSING PERCENT THRESHOLD')
        os._exit(1)
#     elif args.output == None:
#         print('MISSING OUTPUT DIRECTORY')
#         os._exit(1)
    
    # Read in the data
    latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id = ReadATL03(args.input, args.laser)
    
    # Find the epsg code
    epsg_code = convert_wgs_to_utm(longitude[0], latitude[0])
    
    # Orthometrically correct the data using the epsg code
    lat_utm, lon_utm, photon_h = OrthometricCorrection(latitude, longitude, photon_h, epsg_code)
    
    # Get ref-elev and ref_azimuth at photon level
    # Get length of photon array
    heights_len = len(photon_h)
    # Assign segment id to each photon for the segment it is in
    Ph_segment_id = getAtl03SegID(ph_index_beg, segment_id, heights_len)
    # Cast as an int
    Ph_segment_id = Ph_segment_id.astype(np.int)
    # Ref_elev on a per photon level (assign seg ref_elev to photons)
    Ph_ref_elev = ref_elev[np.searchsorted(segment_id, Ph_segment_id)]
    # Ref_azimuth on a per photon level (assign seg ref_azimuth to photons)
    Ph_ref_azimuth = ref_azimuth[np.searchsorted(segment_id, Ph_segment_id)]

    # Aggregate data into dataframe
    dataset_sea = pd.DataFrame({'latitude': lat_utm, 'longitude': lon_utm, 'photon_height': photon_h, 'confidence':conf, 'ref_elevation':Ph_ref_elev, 'ref_azminuth':Ph_ref_azimuth}, 
                           columns=['latitude', 'longitude', 'photon_height', 'confidence', 'ref_elevation', 'ref_azminuth'])

    # Filter data that should not be analyzed
    # Filter for quality flags
    dataset_sea1 = dataset_sea[(dataset_sea.confidence != 0)  & (dataset_sea.confidence != 1)]
    # Filter for elevation range
    dataset_sea1 = dataset_sea1[(dataset_sea1['photon_height'] > -40) & (dataset_sea1['photon_height'] < 5)]
    
    # Focus on specific latitude
    if args.start_lat is not None:
        dataset_sea1 = dataset_sea1[(dataset_sea1['latitude'] > args.start_lat) & (dataset_sea1['latitude'] < args.end_lat)]

    # Bin dataset
    binned_data_sea = bin_data(dataset_sea1, args.lat_res, args.h_res)
    
    # Find mean sea height
    sea_height = get_sea_height(binned_data_sea)

    # Set sea height
    WSHeight = np.nanmedian(sea_height)

    # Calculate sea temperature
    try:
        sst_calculated = get_water_temp(args.input, latitude, longitude)
    
    except OSError:
        sst_calculated = 20
        
    if 15 <= sst_calculated <= 30:
        waterTemp = sst_calculated
    else:
        waterTemp = 20

    # Correct for refraction 
    RefX, RefY, RefZ, RefConf, rawX, rawY, rawZ, ph_ref_azi, ph_ref_elev = RefractionCorrection(waterTemp, WSHeight, 532, dataset_sea1.ref_elevation, dataset_sea1.ref_azminuth, 
                                                                                                dataset_sea1.photon_height, dataset_sea1.longitude, dataset_sea1.latitude, dataset_sea1.confidence)

    # Find bathy depth
    depth = WSHeight - RefZ

    # Create new dataframe with refraction corrected data
    dataset_bath = pd.DataFrame({'latitude': rawY, 'longitude': rawX, 'cor_latitude':RefY, 'cor_longitude':RefX, 'cor_photon_height':RefZ, 'photon_height': rawZ, 'confidence':RefConf, 'depth':depth}, 
                       columns=['latitude', 'longitude', 'photon_height', 'cor_latitude','cor_longitude', 'cor_photon_height', 'confidence', 'depth'])

    # Export dataframe to gpkg
    geodf = geopandas.GeoDataFrame(dataset_bath, geometry=geopandas.points_from_xy(dataset_bath.longitude, dataset_bath.latitude))
    file = args.input[-39:-3]
    geodf.to_file(file + ".gpkg", driver="GPKG")
    
    # Bin dataset again for bathymetry
    binned_data = bin_data(dataset_bath, args.lat_res, args.h_res)
    
    # Find bathymetry 
    bath_height = get_bath_height(binned_data, args.thresh, WSHeight, args.h_res)
    
    # Create figure
    plt.close()
    produce_figures(binned_data, bath_height, sea_height, 10, -20, args.thresh, file)

if __name__ == '__main__':
    main()