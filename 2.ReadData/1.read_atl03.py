import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import pyproj
from pyproj import Proj
from pyproj import Transformer
import argparse


def OrthometricCorrection(lat, lon, Z, epsg):
    # transform ellipsod (WGS84) height to orthometric height
    transformerh = Transformer.from_crs("epsg:4326", "epsg:3855")
    Y_egm08, X_egm08, Z_egm08 = transformerh.transform(lat, lon, Z)
    
    # transform WGS84 proj to local UTM
    myProj = Proj(epsg)
    Y_utm, X_utm = myProj(latitude, longitude)
    
    return Y_utm, X_utm, Z_egm08

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
    
    #ref_elev = f['/' + laser + '/geolocation/ref_elev'][...,]
    #ref_azimuth = f['/' + laser + '/geolocation/ref_azimuth'][...,]
    #ph_index_beg = f['/' + laser + '/geolocation/ph_index_beg'][...,]
    #segment_id = f['/' + laser + '/geolocation/segment_id'][...,]
    
    return latitude, longitude, photon_h

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Specify the input ICESAT H5 file")
    parser.add_argument("-l", "--laser", type=str, help="Specify the ICESAT-2 laser number (1, 2 or 3)")
    parser.add_argument("-e", "--epsg_num", type=int, help="Specify the UTM Zone EPSG code (www.spatialreference.org)")
    #parser.add_argument("-wt", "--waterTemp", type=float, help="Specify the water temperature in degrees C")
    
    args = parser.parse_args()
    
    if args.input == None:
        print('MISSING H5')
        os._exit(1)
    elif args.laser == None:
        print('MISSING LASER NUMBER')
        os._exit(1)
    elif args.epsg_num == None:
        print('MISSING UTM ZONE')
        os._exit(1)

    latitude, longitude, photon_h = ReadATL03(args.input, args.laser)
    
    lat_utm, lon_utm, h_utm = OrthometricCorrection(latitude, longitude, photon_h args.epsg)
    
if __name__ == '__main__':
    main()
