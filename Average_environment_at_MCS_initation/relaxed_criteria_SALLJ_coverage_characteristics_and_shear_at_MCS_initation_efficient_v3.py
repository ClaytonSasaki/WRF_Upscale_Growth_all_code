#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sep 18 2023

@author: crs326

Gets MCS lat/lon centroid and times at initiation for MCSs in a defined region near SDC. 

Finds average environmental characteristics (shear, SALLJ) in a region around MCS initiation location centroid. 

If percentage area coverage of the SALLJ meets a criteria, a SALLJ is found and the average SALLJ characteristics are used.

based on 'relaxed_criteria_SALLJ_spatial_areas_coverage_efficient_v3.py'

"""

import matplotlib
matplotlib.use('Agg') 

from datetime import datetime
import csv
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.dates as mdates
import numpy as np
import numpy.ma as ma
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import pandas as pd
import os
from scipy import stats
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
import math
import xarray

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, vinterp, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

# read in MCS tracks
MCS_tracks_file_path = '/home/disk/monsoon/relampago/analysis//mcs_tracks/wrf/stats/mcs_tracks_pf_20181015_20190430.nc'
MCS_tracks_ncfile = Dataset(MCS_tracks_file_path,'r')

## get MCS initation times and centroid locations
#datetime_str = getvar(MCS_tracks_ncfile, 'datetimestring', meta=False)
#MCS_center_lons = getvar(MCS_tracks_ncfile, 'pf_lon', meta=False)
#MCS_center_lats = getvar(MCS_tracks_ncfile, 'pf_lat', meta=False)

# gets MCS times in byte format
MCS_datetime_bytes = MCS_tracks_ncfile.variables['datetimestring']

# gets MCS iniation times still in byte format
MCS_datetime_initiation_bytes = MCS_datetime_bytes[:,0,:]

# concatenates the bytes that make up one MCS initation time into a single byte string and decodes to a string
MCS_datetime_initiation_str = np.array([bytes(bytes_for_one_time).decode('UTF-8') for bytes_for_one_time in MCS_datetime_initiation_bytes])

# another way to do the same as the line above with the join function instead
#print([b''.join(row).decode('UTF-8') for row in datetime_bytes])

# gets the MCS center lons and lats (nmaxpf chosen to be 0, NOTE: not sure what nmaxpf is?????????)
MCS_center_lons = MCS_tracks_ncfile.variables['pf_lon'][:,:,0]
MCS_center_lats = MCS_tracks_ncfile.variables['pf_lat'][:,:,0]

# get MCS initation center lons and lats
MCS_center_lons_initiation = MCS_center_lons[:,0]
MCS_center_lats_initiation = MCS_center_lats[:,0]

#print(MCS_time_initation_str)
#print(MCS_center_lons_initation)
#print(MCS_center_lats_initation)
print('dates of MCS init \n', MCS_datetime_initiation_str)
print('lons of MCS init \n', MCS_center_lons_initiation)
print('lats of MCS init \n', MCS_center_lats_initiation)
# close the MCS tracks file
MCS_tracks_ncfile.close()

# save MCS initation times and centroid locations of those within defined box near the SDC
lon_min = -66.25
lon_max = -61.1
lat_min = -34.0
lat_max = -29.5

filtered_indices = []

for index, (lat, lon) in enumerate(zip( MCS_center_lats_initiation, MCS_center_lons_initiation)):
    #print('lat', lat)
    #print('lon', lon)
    #print(lat_min <= lat <= lat_max)
    #print(lon_min <= lon <= lon_max)
    if lat_min <= lat <= lat_max and lon_min <= lon <= lon_max:
        #print('yes')
        filtered_indices.append(index)
    else:
        pass
        #print('no')
                                   
MCS_datetime_initiation_filtered_str = MCS_datetime_initiation_str[filtered_indices]
MCS_center_lons_initiation_filtered = MCS_center_lons_initiation[filtered_indices]
MCS_center_lats_initiation_filtered = MCS_center_lats_initiation[filtered_indices]

# number of MCSs within chosen region
num_MCS_init_region = len(MCS_datetime_initiation_filtered_str)
print('# of MCSs within chosen region: ', num_MCS_init_region)

# go through list of times/centroids
for MCS_datetime, MCS_center_lon, MCS_center_lat in zip(MCS_datetime_initiation_filtered_str, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered):

    # get lats/lons of region based on centroid
    lat_bottom_left = MCS_center_lat - 0.75 
    lon_bottom_left = MCS_center_lon - 0.75
                                   
    lat_top_right = MCS_center_lat + 0.75 
    lon_top_right = MCS_center_lon + 0.75

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
    
    MCS_datetime_str = MCS_datetime_dt.strftime('%Y%m%d')

    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(MCS_datetime_dt.strftime('%Y'), MCS_datetime_dt.strftime('%m'), MCS_datetime_dt.strftime('%d'), MCS_datetime_dt.strftime('%H'))

    wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name
    
    print('paths of MCS within region \n', wrf_file_path)

    # get the netCDF
    wrf_ncfile = Dataset(wrf_file_path,'r')
    
    # find if SALLJ is present and if so calculate characteristics
    
    # calculate shear characteristics
    
    # calculate Richardson number
    
    # save MCS characteristics (time, centroid, growth rate, cloud sheild area) and environmental characteristics (SALLJ characterisitcs, shear) to netCDF or dateframe file
    
    wrf_ncfile.close()

    