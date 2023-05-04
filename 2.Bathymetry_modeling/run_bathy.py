#!/usr/bin/env python
# coding: utf-8

'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''


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
from datetime import datetime
# Import functions from utils.py
from bathy_utils import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Specify the input ICESAT H5 file")
    parser.add_argument("-l", "--laser", type=str, help="Specify the ICESAT-2 laser number (1, 2 or 3)")
    parser.add_argument("-th", "--thresh", type=int, help="Specify the threshold percentage")
    parser.add_argument("-o", "--output", type=str, required = False, help="Specify the output location")
    parser.add_argument("-lr", "--lat_res", type=float, default = 10, help="Specify the latitudinal resoltuion (normally 10)")
    parser.add_argument("-hr", "--h_res", type=float, default = 0.5, help="Specify the height resolution (normally 0.5)")
    parser.add_argument("-wt", "--waterTemp", type=float, default = None, required = False, help="Specify the water temperature in degrees C")
    parser.add_argument("-slat", "--start_lat", type=float, required = False, help="Specify the start latitude")
    parser.add_argument("-elat", "--end_lat", type=float, required = False, help="Specify the stop latitude")
    
    args = parser.parse_args()
    
    if args.input == None:
        print('MISSING H5')
        os._exit(1)
    elif args.laser == None:
        print('MISSING LASER NUMBER')
        os._exit(1)
    #elif args.waterTemp == None:
    #    print('MISSING WATER TEMP')
    #    os._exit(1)
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
    
    print('''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.''')
    
    # Read in the data
    latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, alt_sc, seg_ph_count = ReadATL03(args.input, args.laser)
    
    # Find the epsg code
    epsg_code = convert_wgs_to_utm(latitude[0], longitude[0])
    epsg_num = int(epsg_code.split(':')[-1])
    # Orthometrically correct the data using the epsg code
    lat_utm, lon_utm, photon_h = OrthometricCorrection(latitude, longitude, photon_h, epsg_code)
    # count number of photons in each segment: DEPRECATED
    #Ph_num_per_seg = count_ph_per_seg(ph_index_beg, photon_h)
    Ph_num_per_seg = seg_ph_count[ph_index_beg>0]
    # Cast as an int
    Ph_num_per_seg = Ph_num_per_seg.astype(np.int64)
        
    # count_ph_per_seg() function removes zeros from ph_index_beg
    # These 0s are nodata vals in other params (ref_elev etc)
    # Thus no pre-processing is needed as it will map correctly given the nodata values are eliminated
    Ph_ref_elev = ref_linear_interp(Ph_num_per_seg, ref_elev[ph_index_beg>0])
    Ph_ref_azimuth = ref_linear_interp(Ph_num_per_seg, ref_azimuth[ph_index_beg>0])
    Ph_sat_alt = ref_linear_interp(Ph_num_per_seg, alt_sc[ph_index_beg>0]) 
    
    #plt.scatter(np.arange(1, len(Ph_ref_azimuth)+1, 1), Ph_ref_azimuth)
    #plt.show()
    
    
    # Aggregate data into dataframe
    dataset_sea = pd.DataFrame({'latitude': lat_utm, 'longitude': lon_utm, 'photon_height': photon_h, 'confidence':conf, 'ref_elevation':Ph_ref_elev, 'ref_azminuth':Ph_ref_azimuth, 'ref_sat_alt':Ph_sat_alt}, 
                           columns=['latitude', 'longitude', 'photon_height', 'confidence', 'ref_elevation', 'ref_azminuth', 'ref_sat_alt'])
    
    #plt.scatter(lat_utm, photon_h, c='black', s=0.1, alpha=0.1)
    #plt.show()
    # Filter data that should not be analyzed
    # Filter for quality flags
    print('filter quality flags')
    dataset_sea1 = dataset_sea[(dataset_sea.confidence != 0)  & (dataset_sea.confidence != 1)]
    # Filter for elevation range
    dataset_sea1 = dataset_sea1[(dataset_sea1['photon_height'] > -40) & (dataset_sea1['photon_height'] < 5)]
    
    # Focus on specific latitude
    if args.start_lat is not None:
        dataset_sea1 = dataset_sea1[(dataset_sea1['latitude'] > args.start_lat) & (dataset_sea1['latitude'] < args.end_lat)]

    #plt.scatter(dataset_sea1['latitude'], dataset_sea1['photon_height'],c='black',s=0.2,alpha=0.1)
    #plt.show()
    # Bin dataset
    print(dataset_sea1.head())
    binned_data_sea = bin_data(dataset_sea1, args.lat_res, args.h_res)
    
    # Find mean sea height
    sea_height = get_sea_height(binned_data_sea)
    
    # Set sea height
    WSHeight = np.nanmedian(sea_height)

    waterTemp = get_water_temp(args.input, latitude, longitude)
    
    print(waterTemp)
    '''
    # Calculate sea temperature
    if args.waterTemp is not None:
        waterTemp = args.waterTemp
    else:
        try:
            waterTemp = get_water_temp(args.input, latitude, longitude)
        except Exception as e:
            print('NO SST PROVIDED OR RETRIEVED: 20 degrees C assigned')
            waterTemp = 20 
        
    print("water temp:", waterTemp)
    
    #except OSError:
    #    sst_calculated = 20
        
    #if 15 <= sst_calculated <= 30:
    #    waterTemp = sst_calculated
    #else:
    #    waterTemp = 20

    # Correct for refraction 
    print('refrac correction')
    RefX, RefY, RefZ, RefConf, rawX, rawY, rawZ, ph_ref_azi, ph_ref_elev = RefractionCorrection(waterTemp, WSHeight, 532, dataset_sea1.ref_elevation, dataset_sea1.ref_azminuth, dataset_sea1.photon_height, dataset_sea1.longitude, dataset_sea1.latitude, dataset_sea1.confidence, dataset_sea1.ref_sat_alt)
    
    # Find bathy depth
    depth = WSHeight - RefZ

    # Create new dataframe with refraction corrected data
    dataset_bath = pd.DataFrame({'latitude': rawY, 'longitude': rawX, 'cor_latitude':RefY, 'cor_longitude':RefX, 'cor_photon_height':RefZ, 'photon_height': rawZ, 'confidence':RefConf, 'depth':depth}, 
                       columns=['latitude', 'longitude', 'photon_height', 'cor_latitude','cor_longitude', 'cor_photon_height', 'confidence', 'depth'])

     # Export dataframe to gpkg
    #geodf = geopandas.GeoDataFrame(dataset_bath, geometry=geopandas.points_from_xy(dataset_bath.longitude, dataset_bath.latitude))
    file = args.input
    #geodf.to_file(file + ".gpkg", driver="GPKG")
    
    # Bin dataset again for bathymetry
    binned_data = bin_data(dataset_bath, args.lat_res, args.h_res)
    
    # Find bathymetry 
    bath_height, geo_df = get_bath_height(binned_data, args.thresh, WSHeight, args.h_res)
    
    # Create figure
    plt.close()
    produce_figures(binned_data, bath_height, sea_height, 10, -20, args.thresh, file, geo_df, RefY, RefZ, args.laser, epsg_num)
    '''
if __name__ == '__main__':
    main()

