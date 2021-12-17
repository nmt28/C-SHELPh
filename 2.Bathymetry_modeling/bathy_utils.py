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

def ReadATL03(h5_file, laser_num):
    # Read File
    f = h5.File(h5_file,'r')
    
    # Select a laser
    orientation = f['/orbit_info/sc_orient'][0]
    
    # selects the strong beams only [we can include weak beams later on]
    orientDict = {0:'l', 1:'r', 21:'error'}
    laser = 'gt' + laser_num + orientDict[orientation]
    
    # Read in the required photon level data
    photon_h = f['/' + laser + '/heights/h_ph'][...,]
    latitude = f['/' + laser + '/heights/lat_ph'][...,]
    longitude = f['/' + laser + '/heights/lon_ph'][...,]
    conf = f['/' + laser + '/heights/signal_conf_ph/'][...,0]

    # params needed for refraction correction
    
    ref_elev = f['/' + laser + '/geolocation/ref_elev'][...,]
    ref_azimuth = f['/' + laser + '/geolocation/ref_azimuth'][...,]
    ph_index_beg = f['/' + laser + '/geolocation/ph_index_beg'][...,]
    segment_id = f['/' + laser + '/geolocation/segment_id'][...,]
    
    return latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id

# convert_wgs_to_utm function, see https://stackoverflow.com/a/40140326/4556479
def convert_wgs_to_utm(lon: float, lat: float):
    """Based on lat and lng, return best utm epsg-code"""
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = 'epsg:326' + utm_band
        return epsg_code
    epsg_code = 'epsg:327' + utm_band
    return epsg_code

def OrthometricCorrection(lat, lon, Z, epsg):
    # transform ellipsod (WGS84) height to orthometric height
    transformerh = Transformer.from_crs("epsg:4326", "epsg:3855", always_xy=True)
    X_egm08, Y_egm08, Z_egm08 = transformerh.transform(lon, lat, Z)
    
    # transform WGS84 proj to local UTM
    myProj = Proj(epsg)
    X_utm, Y_utm = myProj(lon, lat)
    
    return Y_utm, X_utm, Z_egm08

# Snippet by Eric Guenther (via Amy N.) for assigning photons to a segment
def getAtl03SegID(atl03_ph_index_beg, atl03_segment_id, atl03_heights_len): 
    
    # We need to know spacecraft orbin info, which is provided across segments. This first function assigns photons to the segmenet they belong to. We end up making a new array
    # that has more points to match the photon. Segment is defined as every 100m in the long track.
    
    # Filter all data where atl03_ph_index starts at 0 (0 indicates errors)
    indsNotZero = atl03_ph_index_beg != 0
    atl03_ph_index_beg = atl03_ph_index_beg[indsNotZero]
    atl03_segment_id = atl03_segment_id[indsNotZero]
    
    # Subtract 1 from ph_index_beg to start at python 0th pos
    atl03_ph_index_beg = atl03_ph_index_beg - 1
    
    # Sometimes the ph_index_beg is not at the 0th position, it is is not,
    # add it in and then add the associated segment id
    # Warning, this is assuming that the segment id for the points are from
    # the segment id directly before it, this assumption might fail but I have
    # not come across a case yet where it does.  If you want to play it safe
    # you could comment this section out and then if the first position is not
    # 0 then all photons before the first position will not be assigned a
    # segment id.
    # if atl03_ph_index_beg[0] != 0:
    #     atl03_ph_index_beg = np.append(0,atl03_ph_index_beg)
    #     first_seg_id = atl03_segment_id[0] -1
    #     atl03_segment_id = np.append(first_seg_id,atl03_segment_id)
    
    
    # Append atl03_height_len to end of array for final position
    atl03_ph_index_beg = np.append(atl03_ph_index_beg,atl03_heights_len)
    
    # Make array equal to the length of the atl03_heights photon level data
    ph_segment_id = np.zeros(atl03_heights_len)
    
    # Iterate through ph_index_beg, from the first to second to last number
    # and set the photons between ph_index_beg i to ph_index_beg i + 1 to
    # segment id i
    for i in range(0,len(atl03_ph_index_beg) - 1):
        ph_segment_id[atl03_ph_index_beg[i]:atl03_ph_index_beg[i+1]] = atl03_segment_id[i]
    
    # Return list of segment_id at the photon level
    return ph_segment_id

def bin_data(dataset, lat_res, height_res):
    '''Bin data along vertical and horizontal scales for later segmentation'''
    
    # Calculate number of bins required both vertically and horizontally with resolution size
    lat_bin_number = round(abs(dataset['latitude'].min() - dataset['latitude'].max())/lat_res)
    height_bin_number = round(abs(dataset['photon_height'].min() - dataset['photon_height'].max())/height_res)
    
     # Duplicate dataframe
    dataset1 = dataset
    
    # Cut lat bins
    lat_bins = pd.cut(dataset['latitude'], lat_bin_number, labels = np.array(range(lat_bin_number)))
    
    # Add bins to dataframe
    dataset1['lat_bins'] = lat_bins
    
    # Cut height bins
    height_bins = pd.cut(dataset['photon_height'], height_bin_number, labels = np.round(np.linspace(dataset['photon_height'].min(), dataset['photon_height'].max(), num=height_bin_number), decimals = 1))
    
    # Add height bins to dataframe
    dataset1['height_bins'] = height_bins
    dataset1 = dataset1.reset_index(drop=True)

    return dataset1

def get_sea_height(binned_data):
    '''Calculate mean sea height for easier calculation of depth and cleaner figures'''
    
    # Create sea height list
    sea_height = []
    
    # Group data by latitude
    binned_data_sea = binned_data[(binned_data['photon_height'] > -0.5)] # Filter out subsurface data
    grouped_data = binned_data_sea.groupby(['lat_bins'], group_keys=True)
    data_groups = dict(list(grouped_data))
    
    # Loop through groups and return average sea height
    for k,v in data_groups.items():
        # Create new dataframe based on occurance of photons per height bin
        new_df = pd.DataFrame(v.groupby('height_bins').count())
        
        # Return the bin with the highest count
        largest_h_bin = new_df['latitude'].argmax()
        
        # Select the index of the bin with the highest count
        largest_h = new_df.index[largest_h_bin]
        
        # Calculate the median value of all values within this bin
        lat_bin_sea_median = v.loc[v['height_bins']==largest_h, 'photon_height'].median()
        
        # Append to sea height list
        sea_height.append(lat_bin_sea_median)
        del new_df
        
    # Filter out sea height bin values outside 2 SD of mean.
    mean = np.nanmean(sea_height, axis=0)
    sd = np.nanstd(sea_height, axis=0)
    sea_height_1 = np.where((sea_height > (mean + 2*sd)) | (sea_height < (mean - 2*sd)), np.nan, sea_height).tolist()
    
    return sea_height_1

def get_water_temp(data_path, latitude, longitude):
    '''
    Pull down surface water temperature along the track from the JPL GHRSST opendap website.
    
    The GHRSST data are gridded tiles with dimension 17998 x 35999. 
    To get the specific grid tile of the SST, you must convert from lat, lon coordinates
    to the gridded tile ratio of the SST data product using the coordinates of the IS2 data.
    '''
    # Get date from data filename
    date = data_path[-33:-25]
    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    day_of_year = str(datetime.strptime(date, '%Y%m%d').timetuple().tm_yday)
    # Add zero in front of day of year string
    zero_day_of_year = day_of_year.zfill(3)

    # Calculate ratio of latitude from mid-point of IS2 track
    old_lat = latitude.mean()
    old_lat_min = -90
    old_lat_max = 90
    new_lat_min = 0
    new_lat_max = 17998

    new_lat = round(((old_lat - old_lat_min) / (old_lat_max - old_lat_min)) * 
                    (new_lat_max - new_lat_min) + new_lat_min)

    # Calculate ratio of longitude from mid-point of IS2 track
    old_lon = longitude.mean()
    old_lon_min = -180
    old_lon_max = 180
    new_lon_min = 0
    new_lon_max = 35999

    new_lon = round(((old_lon - old_lon_min) / (old_lon_max - old_lon_min)) * 
                    (new_lon_max - new_lon_min) + new_lon_min)

    # Access the SST data using the JPL OpenDap interface
    url = 'https://opendap.jpl.nasa.gov/opendap/OceanTemperature/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/'\
        + str(year) + '/' + str(zero_day_of_year) + '/' + str(date) \
        + '090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc'
    
    dataset = netCDF4.Dataset(url)
    
    # Access the data and convert the temperature from K to C
    water_temp = dataset['analysed_sst'][0,new_lat, new_lon] - 273.15
    return water_temp

def RefractionCorrection(WTemp, WSmodel, Wavelength, Photon_ref_elev, Ph_ref_azimuth, PhotonZ, PhotonX, PhotonY, Ph_Conf):
    
    '''
    WTemp; there is python library that pulls water temp data 
    WSmodel is the value surface height
    Wavelength is fixed
    '''
    
    # Only process photons below water surface model
    PhotonX = PhotonX[PhotonZ<=WSmodel]
    PhotonY = PhotonY[PhotonZ<=WSmodel]
    Photon_ref_elev = Photon_ref_elev[PhotonZ<=WSmodel]
    Ph_ref_azimuth = Ph_ref_azimuth[PhotonZ<=WSmodel]
    Ph_Conf = Ph_Conf[PhotonZ<=WSmodel]
    PhotonZ = PhotonZ[PhotonZ<=WSmodel]
    
    # water temp for refraction correction
    WaterTemp= WTemp
    
    # Refraction coefficient #
    a = -0.000001501562500
    b = 0.000000107084865
    c = -0.000042759374989
    d = -0.000160475520686
    e = 1.398067112092424
    wl = Wavelength
    
    # refractive index of air
    n1 = 1.00029
    
    # refractive index of water
    n2 = (a*WaterTemp**2) + (b*wl**2) + (c*WaterTemp) + (d*wl) + e
    
    # assumption is 0.25416
    # This example is refractionCoef = 0.25449
    # 1.00029 is refraction of air constant
    CorrectionCoef = (1-(n1/n2))
    #########################
    
    #read photon ref_elev to get theta1
    theta1 = np.pi/2 - Photon_ref_elev
    
    # eq 1. Theta2
    theta2 = np.arcsin(((n1*np.sin(theta1))/n2))
    
    # eq 3. S
    # Approximate water Surface = 1.5
    # D  = raw uncorrected depth
    D = WSmodel - PhotonZ
    
    # For Triangle DTS
    S = D/np.cos(theta1)
    
    # eq 2. R
    R = (S*n1)/n2
    Gamma = (np.pi/2)-theta1
    
    # For triangle RPS
    # phi is an angle needed
    phi = theta1-theta2
    
    # P is the difference between raw and corrected YZ location
    P = np.sqrt(R**2 + S**2 - 2*R*S*np.cos(phi))
    
    # alpha is an angle needed
    alpha = np.arcsin((R*np.sin(phi))/P)
    
    # Beta angle needed for Delta Y an d Delta Z
    Beta = Gamma - alpha
    
    # Delta Y
    DY = P*np.cos(Beta)
    
    # Delta Z
    DZ = P*np.sin(Beta)
    
    # Delta Easting
    DE = DY*np.sin(Ph_ref_azimuth)
    
    # Delta Northing
    DN = DY*np.cos(Ph_ref_azimuth)
    
    outX = PhotonX + DE
    outY = PhotonY + DN
    outZ = PhotonZ + DZ
    '''
        print('\nFor selected Bathy photon:')
        print('lat = ', PhotonY[9000])
        print('long = ', PhotonX[9000])
        print('Raw Depth = ', PhotonZ[9000])
        print('D = ', D[9000])
        
        print('ref_elev = ', Photon_ref_elev[9000])
        
        print('Delta Easting = ', DE[9000])
        print('Delta Northing = ', DN[9000])
        print('Delta Z = ', DZ[9000])
        '''
    return(outX, outY, outZ, Ph_Conf, PhotonX, PhotonY, PhotonZ, Ph_ref_azimuth, Photon_ref_elev) # We are most interested in outx, outy, outz

def get_bath_height(binned_data, percentile, WSHeight, height_resolution):
    '''Calculate bathymetry level per bin based on horizontal resolution'''
    # Create sea height list
    bath_height = []
    
    geo_photon_height = []
    geo_longitude = []
    geo_latitude = []
    
    # Group data by latitude
    # Filter out surface data that are two bins below median surface value calculated above
    binned_data_bath = binned_data[(binned_data['photon_height'] < WSHeight - (height_resolution * 2))] 
    grouped_data = binned_data_bath.groupby(['lat_bins'], group_keys=True)
    data_groups = dict(list(grouped_data))
    
    # Create a percentile threshold of photon counts in each grid, grouped by both x and y axes.
    count_threshold = np.percentile(binned_data.groupby(['lat_bins', 'height_bins']).size().reset_index().groupby('lat_bins')[[0]].max(), percentile)
    
    # Loop through groups and return average bathy height
    for k,v in data_groups.items():
        new_df = pd.DataFrame(v.groupby('height_bins').count())
        bath_bin = new_df['cor_latitude'].argmax()
        bath_bin_h = new_df.index[bath_bin]
        
        # Set threshold of photon counts per bin
        if new_df.iloc[bath_bin]['latitude'] >= count_threshold:
            
            geo_photon_height.append(v.loc[v['height_bins']==bath_bin_h, 'cor_photon_height'].values)
            geo_longitude.append(v.loc[v['height_bins']==bath_bin_h, 'cor_longitude'].values)
            geo_latitude.append(v.loc[v['height_bins']==bath_bin_h, 'cor_latitude'].values)
            
            bath_bin_median = v.loc[v['height_bins']==bath_bin_h, 'cor_photon_height'].median()
            bath_height.append(bath_bin_median)
            del new_df
            
        else:
            bath_height.append(np.nan)
            del new_df
            
    geo_longitude_list = np.concatenate(geo_longitude).ravel().tolist()
    geo_latitude_list = np.concatenate(geo_latitude).ravel().tolist()
    geo_photon_list = np.concatenate(geo_photon_height).ravel().tolist()
    geo_depth = WSHeight - geo_photon_list
    geo_df = pd.DataFrame({'longitude': geo_longitude_list, 'latitude':geo_latitude_list, 'photon_height': geo_photon_list, 'depth':geo_depth})
    
    del geo_longitude_list, geo_latitude_list, geo_photon_list
    
    return bath_height, geo_df

def produce_figures(binned_data, bath_height, sea_height, y_limit_top, y_limit_bottom, percentile, file, geo_df, RefY, RefZ, laser, epsg_num):
    '''Create figures'''
    
    # Create bins for latitude
    x_axis_bins = np.linspace(binned_data.latitude.min(), binned_data.latitude.max(), len(sea_height))
    
    # Create new dataframes for median values  
    bath_median_values = {'x' : x_axis_bins ,'y' : bath_height}
    bath_median_df = pd.DataFrame.from_dict(bath_median_values, orient='index')
    bath_median_df = bath_median_df.transpose()

      
    # Create uniform sea surface based on median sea surface values and filter out surface breaching
    sea_height1 = [np.nanmedian(sea_height) if i == i else np.nan for i in sea_height]
    sea_median_df = pd.DataFrame({'x':x_axis_bins, 'y':sea_height1})

    # Define figure size
    fig = plt.rcParams["figure.figsize"] = (40,5)
    
    # Plot raw points
#     plt.scatter(x=binned_data.latitude, y = binned_data.photon_height, marker='o', lw=0, s=1, alpha = 0.8, c = 'yellow', label = 'Raw photon height')
    plt.scatter(RefY, RefZ, s=0.2, alpha=0.1, c='black')
    plt.scatter(geo_df.latitude, geo_df.photon_height, s=0.5, alpha=0.1, c='red')
    
    #plt.scatter(x=geo_df.latitude, y = geo_df.photon_height, marker='o', lw=0, s=0.8, alpha = 0.8, c = 'black', label = 'Corrected photon bin')

    # Plot median values
    #plt.scatter(bath_median_df.x, bath_median_df.y, marker = 'o', c='r', alpha = 0.8, s = 5, label = 'Median bathymetry')
    plt.scatter(sea_median_df.x, sea_median_df.y, marker = 'o', c='b', alpha = 1, s = 0.5, label = 'Median sea surface')
    
    # Insert titles and sub-titles
    plt.title('Icesat2 Bathymetry\n' + file)
    plt.xlabel('Latitude')
    plt.ylabel('Photon H (in m)')

    plt.legend(loc="upper left")

    # Limit the x and y axes using parameters
    plt.xlim(left=binned_data.latitude.min(), right=binned_data.latitude.max())
    plt.ylim(top = y_limit_top, bottom = y_limit_bottom)
    
    timestr = time.strftime("%Y%m%d_%H%M%S")
    file = file.replace('.h5','')
    # Define where to save file
    plt.savefig(file + '_gt' + str(laser) + '_' + str(percentile) + '_' + timestr + '.png')
    #plt.show()
    #plt.close()
        
    # convert corrected locations back to wgs84 (useful to contain)
    transformer = Transformer.from_crs("EPSG:"+str(epsg_num), "EPSG:4326")

    lon_wgs84, lat_wgs84 = transformer.transform(geo_df.longitude.values, geo_df.longitude.values)

    geo_df['lon_wgs84'] = lon_wgs84
    geo_df['lat_wgs84'] = lat_wgs84
    
    geodf = geopandas.GeoDataFrame(geo_df, geometry=geopandas.points_from_xy(geo_df.lon_wgs84, geo_df.lat_wgs84))
    
    geodf.set_crs(epsg=4326, inplace=True)
    geodf.to_crs(epsg=4326, inplace=True)
    
    geodf.to_file(file + '_gt' + str(laser) + '_' + str(percentile) + ' _EPSG' + str(epsg_num) '_' + timestr + ".gpkg", driver="GPKG")
