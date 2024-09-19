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

from datetime import timedelta
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

# get MCS initation times and centroid locations
datetime_str = getvar(MCS_tracks_ncfile, 'datetimestring')
MCS_center_lons = getvar(MCS_tracks_ncfile, 'pf_lon')
MCS_center_lats = getvar(MCS_tracks_ncfile, 'pf_lat')

MCS_time_initation_str = datetime_str[:,0,0]
MCS_center_lons_initation = MCS_center_lon[:,0,0]
MCS_center_lats_initation = MCS_center_lat[:,0,0]

print(MCS_initation_times_str)

# save MCS initation times and centroid locations of those within defined box near the SDC
lon_min = -66.25
lon_max = -61.1
lat_min = -29.5
lat_max = -34.0

filtered_indices = []

for index, (lat, lon) in enumerate(zip(MCS_time_initation_str, MCS_center_lats_initation, MCS_center_lons_initation):
    if lat_min <= lat <= lat_max and lon_min <= lon <= lon_max:
        filtered_indices.append(index)
                                   
MCS_time_initation_filtered_str = MCS_time_initation_str[filtered_indices]
MCS_center_lons_initation_filtered = MCS_center_lons_initation[filtered_indices]
MCS_center_lats_initation_filtered = MCS_center_lats_initation[filtered_indices]

# go through list of times/centroids
for MCS_time, MCS_center_lon, MCS_center_lat in zip(MCS_time_initation_filtered_str, MCS_center_lons_initation_filtered, MCS_center_lats_initation_filtered):

    # get lats/lons of region based on centroid
    lat_bottom_left = MCS_center_lat - 0.75 
    lon_bottom_left = MCS_center_lon - 0.75
                                   
    lat_top_right = MCS_center_lat + 0.75 
    lon_top_right = MCS_center_lon + 0.75

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_date_str = MCS_time.date().strftime('%Y%m%d')

    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(SALLJ_time.strftime('%Y'), SALLJ_time.strftime('%m'), SALLJ_time.strftime('%d'), SALLJ_time.strftime('%H'))

    wrf_file_path = path_wrf + dt_date_str + wrf_file_name

    # get the netCDF
    wrf_ncfile = Dataset(file_path,'r')
    
    # find if SALLJ is present and if so calculate characteristics
    
    # calculate shear characteristics
    
    # calculate Richardson number
    
    # save MCS characteristics (time, centroid, growth rate, cloud sheild area) and environmental characteristics (SALLJ characterisitcs, shear) to netCDF or dateframe file

    