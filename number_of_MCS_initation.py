#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sep 18 2023

@author: crs326

Gets MCS lat/lon centroid and times at initiation for 'new track' MCSs in a defined region near SDC. 

based on 'relaxed_criteria_SALLJ_coverage_characteristics_and_shear_at_MCS_initation_efficient_v3.py'

"""

import matplotlib
matplotlib.use('Agg') 

from datetime import datetime, timedelta
import csv
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.dates as mdates
import numpy as np
import numpy.ma as ma
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import pickle
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

offset_MCS_and_conditions = True

hours_offset = -1

######### get MCS initation times and centroid locations #########

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

MCS_start_type_location_filtered = np.array(MCS_tracks_ncfile.variables['starttrackresult'])

print(MCS_start_type_location_filtered)

MCS_start_type_filtered_indices = np.where(MCS_start_type_location_filtered == 10)[0]

MCS_datetime_init_filt_by_loc_and_start_type_str = MCS_datetime_initiation_str[MCS_start_type_filtered_indices]
MCS_center_lons_init_filt_by_loc_and_start_type = MCS_center_lons_initiation[MCS_start_type_filtered_indices]
MCS_center_lats_init_filt_by_loc_and_start_type = MCS_center_lats_initiation[MCS_start_type_filtered_indices]

num_MCS_init_region_and_new_MCS = len(MCS_datetime_init_filt_by_loc_and_start_type_str)
print('# of MCSs within chosen region and new MCS: ', num_MCS_init_region_and_new_MCS)

##### compare MCS times to SALLJ times ###
#
#print('dates of MCS init filtered by location and start type: ', MCS_datetime_init_filt_by_loc_and_start_type_str)
#
#SALLJ_time_list = pickle.load(open("/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/Composite_shear_maps_by_SALLJ_presence/3relaxed_SALLJ_times_%s_%s_%s.dat" %('VMRS', '20181101', '20190430'), "rb"))
#
#SALLJ_time_list_arr = np.array([str(SALLJ_time) for SALLJ_time in SALLJ_time_list])
#
#SALLJ_time_list_arr_dt = [datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in SALLJ_time_list_arr]
#
#print('dates of SALLJ at VMRS: ', SALLJ_time_list_arr_dt)

# close the MCS tracks file
MCS_tracks_ncfile.close()