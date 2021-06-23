import numpy as np
import h5py as h5
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import matplotlib
matplotlib.use("TkAgg")
import pyproj
from pyproj import Proj
from pyproj import Transformer
import pandas as pd
import argparse
import os
import subprocess

# Snippet by Eric Guenther (via Amy N.) for assigning photons to a segment
def getAtl03SegID(atl03_ph_index_beg, atl03_segment_id, atl03_heights_len):
    
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
    
def RefractionCorrection(WTemp, WSmodel, Wavelength, Photon_ref_elev, Ph_ref_azimuth, PhotonZ, PhotonX, PhotonY, Ph_Conf):
    
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
    return(outX, outY, outZ, Ph_Conf, PhotonX, PhotonY, PhotonZ, Ph_ref_azimuth, Photon_ref_elev)


# Read in the required photon level data
Photon_h = f['/' + laser + '/heights/h_ph'][...,]
latitude = f['/' + laser + '/heights/lat_ph'][...,]
longitude = f['/' + laser + '/heights/lon_ph'][...,]
conf = f['/' + laser + '/heights/signal_conf_ph/'][...,0]

ref_elev = f['/' + laser + '/geolocation/ref_elev'][...,]
ref_azimuth = f['/' + laser + '/geolocation/ref_azimuth'][...,]
ph_index_beg = f['/' + laser + '/geolocation/ph_index_beg'][...,]
segment_id = f['/' + laser + '/geolocation/segment_id'][...,]



###### Get ref-elev and ref_azimuth at photon level
# Get length of photon array
heights_len = len(Zegm)
# Assign segment id to each photon for the segment it is in
Ph_segment_id = getAtl03SegID(ph_index_beg, segment_id, heights_len)
# Cast as an int
Ph_segment_id = Ph_segment_id.astype(np.int)
# Ref_elev on a per photon level (assign seg ref_elev to photons)
Ph_ref_elev = ref_elev[np.searchsorted(segment_id, Ph_segment_id)]
# Ref_azimuth on a per photon level (assign seg ref_azimuth to photons)
Ph_ref_azimuth = ref_azimuth[np.searchsorted(segment_id, Ph_segment_id)]

# Input Water surface height until can do it manually
#WSHeight = float(input("What is the WSM Height?"))

# Auto method
WSHeight = WaterSurface

RefX, RefY, RefZ, RefConf, rawX, rawY, rawZ, ph_ref_azi, ph_ref_elev = RefractionCorrection(args.waterTemp, WSHeight, 532, Ph_ref_elev, Ph_ref_azimuth, Zegm, Xutm, Yutm, conf)


