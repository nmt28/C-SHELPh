#!/usr/bin/env python
# coding: utf-8

'''
THE SOFTWARE IS pROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXpRESS OR IMpLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A pARTICULAR pURpOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COpYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''


#####
# This group of functions processes ICESAT2 data and creates a bathymetric model.
# To do this, it follows a number of steps in the form of functions, including:
# 1. Reading data (ReadATL03())
# 2. Orthometrically correcting the dataset (OrthometricCorrection())
# 3. pulling down the data segment ID (getAtl03SegID())
# 4. Bin the data along latitudinal and height gradients (bin_data())
# 5. Calculate sea height (get_sea_height())
# 6. Get water temperature (get_water_temp())
# 7. Correct bathymetry surface for refraction (RefractionCorrection())
# 8. Calculate bathymetry height (get_bath_height())
# 9. produce figures (produce_figures())
#####

import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import pyproj
from pyproj import Proj
from pyproj import Transformer
import pandas as pd
import argparse
import os
import time
import utm
import math
import fiona
import geopandas
from datetime import datetime
# Import functions from utils.py
from bathy_utils import *
# need s3fs installed

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Specify the input ICESAT H5 file")
    parser.add_argument("-l", "--laser", type=str, help="Specify the ICESAT-2 laser number (1, 2 or 3)")
    parser.add_argument("-th", "--thresh", type=int, default=None, help="Specify the threshold percentage")
    parser.add_argument("-tl", "--threshlist", nargs='+', default=None, help="Specify a list of thresholds to trial")
    parser.add_argument("-o", "--output", type=str, required = False, help="Specify the output location")
    parser.add_argument("-lr", "--lat_res", type=float, default = 10, help="Specify the latitudinal resoltuion (normally 10)")
    parser.add_argument("-hr", "--h_res", type=float, default = 0.5, help="Specify the height resolution (normally 0.5)")
    parser.add_argument("-wt", "--water_temp", type=float, default = None, required = False, help="Specify the water temperature in degrees C")
    parser.add_argument("-slat", "--start_lat", type=float, required = False, help="Specify the start latitude")
    parser.add_argument("-elat", "--end_lat", type=float, required = False, help="Specify the stop latitude")
    
    args = parser.parse_args()
    
    if args.input == None:
        print('MISSING H5')
        os._exit(1)
    elif args.laser == None:
        print('MISSING LASER NUMBER')
        os._exit(1)
    #elif args.water_temp == None:
    #    print('MISSING WATER TEMp')
    #    os._exit(1)
    elif args.lat_res == None:
        print('MISSING LATITUDINAL RESOLUTION')
        os._exit(1)
    elif args.h_res == None:
        print('MISSING HEIGHT RESOLUTION')
        os._exit(1)
    elif args.thresh == None:
        if args.threshlist == None:
            print('MISSING pERCENT THRESHOLD VALUE OR LIST')
            os._exit(1)
#     elif args.output == None:
#         print('MISSING OUTpUT DIRECTORY')
#         os._exit(1)
    
    print('''
THE SOFTWARE IS pROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXpRESS OR IMpLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A pARTICULAR pURpOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COpYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.''')
    
    # Read in the data
    latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, alt_sc, seg_ph_count = read_atl03(args.input, args.laser)
    
    # Find the epsg code
    epsg_code = convert_wgs_to_utm(latitude[0], longitude[0])
    epsg_num = int(epsg_code.split(':')[-1])
    # Orthometrically correct the data using the epsg code
    lat_utm, lon_utm, photon_h = orthometric_correction(latitude, longitude, photon_h, epsg_code)
    # count number of photons in each segment: DEpRECATED
    #ph_num_per_seg = count_ph_per_seg(ph_index_beg, photon_h)
    ph_num_per_seg = seg_ph_count[ph_index_beg>0]
    # Cast as an int
    ph_num_per_seg = ph_num_per_seg.astype(np.int64)
        
    # count_ph_per_seg() function removes zeros from ph_index_beg
    # These 0s are nodata vals in other params (ref_elev etc)
    # Thus no pre-processing is needed as it will map correctly given the nodata values are eliminated
    ph_ref_elev = ref_linear_interp(ph_num_per_seg, ref_elev[ph_index_beg>0])
    ph_ref_azimuth = ref_linear_interp(ph_num_per_seg, ref_azimuth[ph_index_beg>0])
    ph_sat_alt = ref_linear_interp(ph_num_per_seg, alt_sc[ph_index_beg>0])

    # Aggregate data into dataframe
    dataset_sea = pd.DataFrame({'latitude': lat_utm, 'longitude': lon_utm, 'photon_height': photon_h, 'confidence':conf, 'ref_elevation':ph_ref_elev, 'ref_azminuth':ph_ref_azimuth, 'ref_sat_alt':ph_sat_alt},
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
    med_water_surface_h = np.nanmedian(sea_height)

    # Calculate sea temperature
    if args.water_temp is not None:
        water_temp = args.water_temp
    else:
        try:
            water_temp = get_water_temp(args.input, latitude, longitude)
        except Exception as e:
            print('NO SST pROVIDED OR RETRIEVED: 20 degrees C assigned')
            water_temp = 20
    
    print("water temp:", water_temp)


    # Correct for refraction 
    print('refrac correction')
    ref_x, ref_y, ref_z, ref_conf, raw_x, raw_y, raw_z, ph_ref_azi, ph_ref_elev = refraction_correction(water_temp, med_water_surface_h, 532, dataset_sea1.ref_elevation, dataset_sea1.ref_azminuth, dataset_sea1.photon_height, dataset_sea1.longitude, dataset_sea1.latitude, dataset_sea1.confidence, dataset_sea1.ref_sat_alt)
    
    # Find bathy depth
    depth = med_water_surface_h - ref_z

    # Create new dataframe with refraction corrected data
    dataset_bath = pd.DataFrame({'latitude': raw_y, 'longitude': raw_x, 'cor_latitude':ref_y, 'cor_longitude':ref_x, 'cor_photon_height':ref_z, 'photon_height': raw_z, 'confidence':ref_conf, 'depth':depth},
                       columns=['latitude', 'longitude', 'photon_height', 'cor_latitude','cor_longitude', 'cor_photon_height', 'confidence', 'depth'])

     # Export dataframe to gpkg
    #geodf = geopandas.GeoDataFrame(dataset_bath, geometry=geopandas.points_from_xy(dataset_bath.longitude, dataset_bath.latitude))
    file = args.input
    #geodf.to_file(file + ".gpkg", driver="GpKG")
    
    # Bin dataset again for bathymetry
    binned_data = bin_data(dataset_bath, args.lat_res, args.h_res)
    
    if isinstance(args.thresh, int) == True:
        # Find bathymetry
        bath_height, geo_df = get_bath_height(binned_data, args.thresh, med_water_surface_h, args.h_res)
        
        # Create figure
        plt.close()
        produce_figures(binned_data, bath_height, sea_height, 10, -20, args.thresh, file, geo_df, ref_y, ref_z, args.laser, epsg_num)
    elif isinstance(args.threshlist, list)==True:
        for thresh in args.threshlist:
            print("using threshold:", str(thresh))
            bath_height, geo_df = get_bath_height(binned_data, int(thresh), med_water_surface_h, args.h_res)
        
            # Create figure
            plt.close()
            produce_figures(binned_data, bath_height, sea_height, 10, -20, str(thresh), file, geo_df, ref_y, ref_z, args.laser, epsg_num)

    
if __name__ == '__main__':
    main()

