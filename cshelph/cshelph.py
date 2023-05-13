#!/usr/bin/env python
# coding: utf-8

'''
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXpRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'''


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
#import netCDF4
from datetime import datetime
import utm
import xarray as xr
import fsspec
import earthaccess
import geopandas

# need s3fs installed

def read_atl03(h5_file, laser_num):
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
    altitude_sc = f['/' + laser + '/geolocation/altitude_sc'][...,]
    seg_ph_count = f['/' + laser + '/geolocation/segment_ph_cnt'][...,]
    
    return latitude, longitude, photon_h, conf, ref_elev, ref_azimuth, ph_index_beg, segment_id, altitude_sc, seg_ph_count
    
def convert_wgs_to_utm(lat, lon):
	easting, northing, num, letter = utm.from_latlon(lat, lon)
	if letter >= 'N':
		epsg = 'epsg:326' + str(num)
	elif letter < 'N':
		epsg = 'epsg:327' + str(num)
	else:
		print('Error Finding UTM')
		
	return epsg

def orthometric_correction(lat, lon, Z, epsg):
    # transform ellipsod (WGS84) height to orthometric height
    transformerh = Transformer.from_crs("epsg:4326", "epsg:3855", always_xy=True)
    X_egm08, Y_egm08, Z_egm08 = transformerh.transform(lon, lat, Z)
    
    # transform WGS84 proj to local UTM
    myproj = Proj(epsg)
    X_utm, Y_utm = myproj(lon, lat)
    
    return Y_utm, X_utm, Z_egm08

    
def count_ph_per_seg(ph_index_beg, photon_h): # DEpRECATED
    
    ph_index_beg = ph_index_beg[ph_index_beg!=0]
    
    # add an extra val at the end of array for upper bounds
    ph_index_beg = np.hstack((ph_index_beg, len(photon_h)+1))-1
    
    photon_id = []
    #iterate over the photon indexes (ph_index_beg)
    for i, num in enumerate(np.arange(0, len(ph_index_beg)-1, 1)):
        photon_id.append(len(photon_h[ph_index_beg[i]:ph_index_beg[i+1]]))
    photon_id = np.array(photon_id)
    
    return photon_id

    
def ref_linear_interp(photon_count, ref_elev):

    arr = []
    for i in range(len(ref_elev)):
        try:
            min = ref_elev[i-1]
            max = ref_elev[i]
        except:
            min = ref_elev[i]
            max = ref_elev[i]
            
        try:
            min = ref_elev[i]
            max = ref_elev[i+1]
        except:
            min = ref_elev[i]
            max = ref_elev[i]
        
        if min==max:
            sub = np.full((photon_count[i]), min)
            arr.append(sub)
        else:
            sub_tmp = np.linspace(min, max, photon_count[i]+1)
            if len(sub_tmp)>1:
                sub = np.linspace(sub_tmp[1], sub_tmp[-1], photon_count[i])
                arr.append(sub)
            else:
                arr.append(sub_tmp)

    return np.concatenate(arr, axis=None).ravel()

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

def get_sea_height(binned_data, surface_buffer):
    '''Calculate mean sea height for easier calculation of depth and cleaner figures'''
    
    # Create sea height list
    sea_height = []
    
    # Group data by latitude
    binned_data_sea = binned_data[(binned_data['photon_height'] > surface_buffer)] # Filter out subsurface data
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

    try:
        auth = earthaccess.login(strategy="netrc")
    except:
        print("Login credentials not found. Sign in manually, your netrc file will be created for next time")
        auth = earthaccess.login(strategy="interactive", persist=True)
    
    # Get date from data filename
    file_date = data_path[-33:-25]
    
    date_range = file_date[0:4] + '-' + file_date[4:6] + '-' + file_date[6:]
    

    location_df = pd.DataFrame({'longitude':longitude,'latitude':latitude})
    
    location_df = location_df.dropna(axis=0)
    
    med_lon = np.nanmedian(location_df['longitude'])
    med_lat = np.nanmedian(location_df['latitude'])
    
    minx, miny, maxx, maxy = list(location_df.total_bounds)
    
    # The first step is to create a DataCollections query

    Query = earthaccess.collection_query()

    # Use chain methods to customize our query
    Query.keyword('GHRSST Level 4 CMC0.1deg Global Foundation Sea Surface Temperature Analysis').bounding_box(minx,miny,maxx,maxy).temporal(date_range,date_range)

    collections = Query.fields(['ShortName','Version']).get(10)

    short_name = collections[0]["umm"]["ShortName"]
    
    Query = earthaccess.granule_query().short_name(short_name).version("3.0").bounding_box(minx,miny,maxx,maxy).temporal("2020-01-01","2020-01-01")
    
    granules = Query.get(10)
    
    ds_L3 = xr.open_mfdataset(
    earthaccess.open(granules),
    combine='nested',
    concat_dim='time',
    coords='minimal',
    )
    
    sea_temp = ds_L3['analysed_sst'].sel(lat=lat_med,lon=lon_med,method='nearest').load()
    
    sst = round(np.nanmedian(sea_temp.values)-273,2)
    
    return sst

def refraction_correction(water_temp, water_surface, wavelength, photon_ref_elev, ph_ref_azimuth, photon_z, photon_x, photon_y, ph_conf, satellite_altitude):
    
    '''
    WTemp; there is python library that pulls water temp data
    water_surface is the value surface height
    Wavelength is fixed
    '''
    
    # Only process photons below water surface model
    photon_x = photon_x[photon_z<=water_surface]
    photon_y = photon_y[photon_z<=water_surface]
    photon_ref_elev = photon_ref_elev[photon_z<=water_surface]
    satellite_altitude = satellite_altitude[photon_z<=water_surface]
    ph_ref_azimuth = ph_ref_azimuth[photon_z<=water_surface]
    ph_conf = ph_conf[photon_z<=water_surface]
    photon_z = photon_z[photon_z<=water_surface]
    
    # Refraction coefficient #
    a = -0.000001501562500
    b = 0.000000107084865
    c = -0.000042759374989
    d = -0.000160475520686
    e = 1.398067112092424
    wl = wavelength
    
    # refractive index of air
    n1 = 1.00029
    
    # refractive index of water
    n2 = (a*water_temp**2) + (b*wl**2) + (c*water_temp) + (d*wl) + e
    
    # assumption is 0.25416
    # This example is refractionCoef = 0.25449
    # 1.00029 is refraction of air constant
    #correction_coef = (1-(n1/n2))
    #########################
    
    # read photon ref_elev to get theta1
    # Does not account for curvature of Earth
    theta1 = np.pi/2 - photon_ref_elev
    
    # H = orbital altitude of IS2 (496km as mean)
    # H = 496km. we pass in the mean of the orbit from /geolocation/altitude_sc/
    # Diff from min to max of 100m over an orbit is 0.02% at 496km
    # More error probably introduced from Re (mean Earth radius) than intra-orbit changes in altitude
    # H = 496
    H = satellite_altitude/1000
    # Re = Radius of Earth (6371km mean)
    Re = 6371
    
    # remove as potentially more inaccurate of a correcttion
    #theta1 = np.arctan((H*np.tan(theta_1))/Re)
    
    # eq 1. Theta2
    theta2 = np.arcsin(((n1*np.sin(theta1))/n2))
    
    # eq 3. S
    # Approximate water Surface = 1.5
    # D  = raw uncorrected depth
    D = water_surface - photon_z
    
    # For Triangle DTS
    S = D/np.cos(theta1)
    
    # eq 2. R
    R = (S*n1)/n2
    Gamma = (np.pi/2)-theta1
    
    # For triangle RpS
    # phi is an angle needed
    phi = theta1-theta2
    
    # p is the difference between raw and corrected YZ location
    p = np.sqrt(R**2 + S**2 - 2*R*S*np.cos(phi))
    
    # alpha is an angle needed
    alpha = np.arcsin((R*np.sin(phi))/p)
    
    # Beta angle needed for Delta Y an d Delta Z
    Beta = Gamma - alpha
    
    # Delta Y
    DY = p*np.cos(Beta)
    
    # Delta Z
    DZ = p*np.sin(Beta)
    
    # Delta Easting
    DE = DY*np.sin(ph_ref_azimuth)
    
    # Delta Northing
    DN = DY*np.cos(ph_ref_azimuth)
    
    out_x = photon_x + DE
    out_y = photon_y + DN
    out_z = photon_z + DZ

    return(out_x, out_y, out_z, ph_conf, photon_x, photon_y, photon_z, ph_ref_azimuth, photon_ref_elev) # We are most interested in out_x, out_y, out_z

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
    x_bins = np.linspace(binned_data.latitude.min(), binned_data.latitude.max(), len(sea_height))
    
    # Create new dataframes for median values  
    bath_median_df = pd.DataFrame.from_dict({'x' : x_bins ,'y' : bath_height}, orient='index')
    bath_median_df = bath_median_df.transpose()

      
    # Create uniform sea surface based on median sea surface values and filter out surface breaching
    sea_surf = [np.nanmedian(sea_height) if i == i else np.nan for i in sea_height]
    sea_median_df = pd.DataFrame({'x':x_bins, 'y':sea_surf})

    # Define figure size
    fig = plt.rcParams["figure.figsize"] = (40,5)
    
    # plot raw points
    plt.scatter(x=binned_data.latitude, y = binned_data.photon_height, marker='o', lw=0, s=1, alpha = 0.8, c = 'yellow', label ='Raw photon height')
    plt.scatter(RefY, RefZ, s=0.2, alpha=0.1, c='black')
    plt.scatter(geo_df.latitude, geo_df.photon_height, s=0.5, alpha=0.1, c='red', label = 'Classified photons')
    
    #plt.scatter(x=geo_df.latitude, y = geo_df.photon_height, marker='o', lw=0, s=0.8, alpha = 0.8, c = 'black', label = 'Corrected photon bin')

    # plot median values
    #plt.scatter(bath_median_df.x, bath_median_df.y, marker = 'o', c='r', alpha = 0.8, s = 5, label = 'Median bathymetry')
    plt.scatter(sea_median_df.x, sea_median_df.y, marker = 'o', c='b', alpha = 1, s = 0.5, label = 'Median sea surface')
    
    # Insert titles and sub-titles
    plt.title('Icesat2 Bathymetry\n' + file)
    plt.xlabel('Latitude', fontsize=25)
    plt.ylabel('photon Height (m)', fontsize=25)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.legend(loc="upper left",prop={'size': 20})

    # Limit the x and y axes using parameters
    plt.xlim(left=binned_data.latitude.min(), right=binned_data.latitude.max())
    plt.ylim(top = y_limit_top, bottom = y_limit_bottom)
    
    timestr = time.strftime("%Y%m%d%H%M%S")
    file = file.replace('.h5','')
    # Define where to save file
    plt.tight_layout()
    plt.savefig(file + '_gt' + str(laser) + '_' + str(percentile) + '_EPSG' + str(epsg_num) + '_' + timestr + ".png")
    #plt.show()
    #plt.close()
        
        # convert corrected locations back to wgs84 (useful to contain)
    transformer = Transformer.from_crs("EPSG:"+str(epsg_num), "EPSG:4326", always_xy=True)
    #print(transformer)
    lon_wgs84, lat_wgs84 = transformer.transform(geo_df.longitude.values, geo_df.latitude.values)

    geo_df['lon_wgs84'] = lon_wgs84
    geo_df['lat_wgs84'] = lat_wgs84
    
    geodf = geopandas.GeoDataFrame(geo_df, geometry=geopandas.points_from_xy(geo_df.lon_wgs84,geo_df.lat_wgs84))
    
    geodf.set_crs(epsg=4326, inplace=True)
    
    geodf.to_file(file + '_gt' + str(laser) + '_' + str(percentile) + '_EPSG' + str(epsg_num) + '_' + timestr + ".gpkg", driver="GPKG")
