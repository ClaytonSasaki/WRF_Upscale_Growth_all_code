#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sep 18 2023

@author: crs326

Gets MCS lat/lon centroid and times at initiation for MCSs in a defined region near SDC. 

Finds average environmental characteristics (shear, SALLJ) in a region around MCS initiation location centroid. 

If percentage area coverage of the SALLJ meets a criteria, a SALLJ is found and the average SALLJ characteristics are used.

Originally based on 'relaxed_criteria_SALLJ_spatial_areas_coverage_efficient_v3.py', 'Find_all_SALLJ_times_relaxed_efficient_avg_SALLJ_level.py', and 'SALLJ_spatial_areas_coverage.py'

ALMOST DUPLICATE OF 'relaxed_criteria_SALLJ_coverage_characteristics_and_shear_at_MCS_initation_efficient_v3.py' but using a different MCS_init_area

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

MCS_init_location_filter = True

MCS_start_type_filter = True

offset_MCS_and_conditions = False

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
MCS_center_lons = MCS_tracks_ncfile.variables['meanlon']
MCS_center_lats = MCS_tracks_ncfile.variables['meanlat']

# get MCS initation center lons and lats
MCS_center_lons_initiation = MCS_center_lons[:,0]
MCS_center_lats_initiation = MCS_center_lats[:,0]

#print(MCS_time_initation_str)
#print(MCS_center_lons_initation)
#print(MCS_center_lats_initation)
#print('dates of MCS init \n', MCS_datetime_initiation_str)
#print('lons of MCS init \n', MCS_center_lons_initiation)
#print('lats of MCS init \n', MCS_center_lats_initiation)

######### create filter MCS tracks by centroid location ######### 
MCS_init_area = 'Zhang_area'

# save MCS initation times and centroid locations of those within defined box near the SDC (ZHANG ET AL. 2021 area)
lon_min = -70.0
lon_max = -59.0
lat_min = -36.0
lat_max = -26.5
        
condition_lat_min = lat_min <= MCS_center_lats_initiation
condition_lat_max = MCS_center_lats_initiation <= lat_max
condition_lon_min = lon_min <= MCS_center_lons_initiation
condition_lon_max = MCS_center_lons_initiation <= lon_max

######### creat filter MCS tracks by MCS start type ######### 

MCS_start_type = np.array(MCS_tracks_ncfile.variables['starttrackresult'])

condition_start_type = MCS_start_type == 10

######### fitler MCS tracks by chosen filters ######### 

if(MCS_init_location_filter == True) and (MCS_start_type_filter == False):

    mask = condition_lat_min & condition_lat_max & condition_lon_min & condition_lon_max
    
    filter_label = '_filtered_init_loc'
    
if(MCS_init_location_filter == False) and (MCS_start_type_filter == True):

    mask = condition_start_type
    
    filter_label = '_filtered_start_type'
    
if(MCS_init_location_filter == True) and (MCS_start_type_filter == True):

    mask = condition_lat_min & condition_lat_max & condition_lon_min & condition_lon_max & condition_start_type
    
    filter_label = '_filtered_init_loc_and_start_type'
    
if(MCS_init_location_filter == False) and (MCS_start_type_filter == False):

    mask = np.full(shape=len(MCS_datetime_initiation_str), fill_value=True, dtype=bool)
    
    filter_label = ''
                                   
MCS_datetime_initiation_str_filtered = MCS_datetime_initiation_str[mask]
MCS_center_lons_initiation_filtered = MCS_center_lons_initiation[mask]
MCS_center_lats_initiation_filtered = MCS_center_lats_initiation[mask]

# number of MCSs that meet chosen conditions
num_MCS = len(MCS_datetime_initiation_str_filtered)
print('# of MCSs that meet chosen conditions: ', num_MCS)

MCS_duration = np.array(MCS_tracks_ncfile.variables['length'])

print('MCS_duration', MCS_duration)
print('len(MCS_duration)', len(MCS_duration))

MCS_duration_filtered = MCS_duration[mask]

MCS_majoraxislength_init = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,0]

MCS_majoraxislength_init_2 = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,1]

MCS_majoraxislength_growth = MCS_majoraxislength_init_2 - MCS_majoraxislength_init

MCS_majoraxislength_growth_filtered = MCS_majoraxislength_growth[mask]

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

######### Get area average environmental data centered on each MCS track chosen ######### 

mean_bulk_shear_0_3km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
mean_bulk_shear_0_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
mean_bulk_shear_2_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
prop_area_SALLJ = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
median_SALLJ_max_wind = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
median_SALLJ_height = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)

SALLJ_count = 0
low_cov_SALLJ_count = 0
med_cov_SALLJ_count = 0
high_cov_SALLJ_count = 0

# go through list of times/centroids for each MCS to get corresponding environmental conditions
for count, (MCS_datetime, MCS_center_lon, MCS_center_lat) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered)):

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
    
    if offset_MCS_and_conditions == True:
        
        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)
        
    else: #offset_MCS_and_conditions == False
        
        conditions_datetime_dt = MCS_datetime_dt
    
    MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')

    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))

    wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name
    
    print('path of MCS within region: ', wrf_file_path)

    # get the netCDF
    wrf_ncfile = Dataset(wrf_file_path,'r')
    
    ####### get environmental profiles #######
    
    # get lats/lons of region based on centroid
    lat_bottom_left_env = MCS_center_lat - 0.75 
    lon_bottom_left_env = MCS_center_lon - 0.75
                                   
    lat_top_right_env = MCS_center_lat + 0.75 
    lon_top_right_env = MCS_center_lon + 0.75
    
    #print(MCS_center_lat, MCS_center_lon)
    #print(lat_bottom_left_env, lon_bottom_left_env)
    #print(lat_top_right_env, lon_top_right_env)
    
    # get xy for regions
    bottom_left_xy = ll_to_xy(wrf_ncfile, lat_bottom_left_env, lon_bottom_left_env)
    top_right_xy = ll_to_xy(wrf_ncfile, lat_top_right_env, lon_top_right_env)  

    # read in the variables 
    hght_original = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in m NOTE: this is AGL not MSL!!
    u_original = getvar(wrf_ncfile, 'ua', units='kt')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts
    v_original = getvar(wrf_ncfile, 'va', units='kt')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts
    
    #print(u_original)
    #print(hght_original)
    
    u2km = interplevel(u_original, hght_original, 2000)
    v2km = interplevel(v_original, hght_original, 2000)

    u3km = interplevel(u_original, hght_original, 3000)
    v3km = interplevel(v_original, hght_original, 3000)
    
    u6km = interplevel(u_original, hght_original, 6000)
    v6km = interplevel(v_original, hght_original, 6000)
    
    u_sfc = u_original[3, :, :]
    v_sfc = v_original[3, :, :]

    # calculate mean values
    mean_u2km = np.nanmean(u2km)
    mean_v2km = np.nanmean(v2km)
    
    mean_u3km = np.nanmean(u3km)
    mean_v3km = np.nanmean(v3km)
    
    mean_u6km = np.nanmean(u6km)
    mean_v6km = np.nanmean(v6km)
    
    mean_u_sfc = np.nanmean(u_sfc)
    mean_v_sfc = np.nanmean(v_sfc)
    
#    mean_u_sfc = u_mean_profile[3]
#    mean_v_sfc = v_mean_profile[3]
#    
#    mean_u2km = interplevel(u_mean_profile, hght_mean_profile, 2000)
#    mean_v2km = interplevel(v_mean_profile, hght_mean_profile, 2000)
#
#    mean_u3km = interplevel(u_mean_profile, hght_mean_profile, 3000)
#    mean_v3km = interplevel(v_mean_profile, hght_mean_profile, 3000)
#    
#    mean_u6km = interplevel(u_mean_profile, hght_mean_profile, 6000)
#    mean_v6km = interplevel(v_mean_profile, hght_mean_profile, 6000)
    
    # calcualte mean shear
    mean_udiff_0_3 = mean_u3km - mean_u_sfc
    mean_vdiff_0_3 = mean_v3km - mean_v_sfc
    
    mean_shear3km = np.sqrt(mean_udiff_0_3**2 + mean_vdiff_0_3**2)

    mean_bulk_shear_0_3km[count] = mean_shear3km
    
    mean_udiff_0_6 = mean_u6km - mean_u_sfc
    mean_vdiff_0_6 = mean_v6km - mean_v_sfc
    
    mean_shear6km = np.sqrt(mean_udiff_0_6**2 + mean_vdiff_0_6**2)

    mean_bulk_shear_0_6km[count] = mean_shear6km
    
    mean_udiff_2_6 = mean_u6km - mean_u2km
    mean_vdiff_2_6 = mean_v6km - mean_v2km
    
    mean_shear2_6km = np.sqrt(mean_udiff_2_6**2 + mean_vdiff_2_6**2)

    mean_bulk_shear_2_6km[count] = mean_shear2_6km
    
    ########## find if SALLJ is present and if so calculate characteristics #########
    
    # for finding LLJs #

    # 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
    crit = [2 ,3200, 5700, 19.4384, 7.77538] # relaxed to 10 m/s max and 4 m/s decrease

    crit_num = crit[0]
    max_search_hgt = crit[1]
    min_search_hgt = crit[2]
    max_wind_threshold = crit[3]
    decrease_to_min_threshold = crit[4]
    
    lat_bottom_left_SALLJ = -30
    lon_bottom_left_SALLJ = -65
    lat_top_right_SALLJ = -28
    lon_top_right_SALLJ = -60
    
    # get xy for regions
    bottom_left_xy = ll_to_xy(wrf_ncfile, lat_bottom_left_SALLJ, lon_bottom_left_SALLJ)
    top_right_xy = ll_to_xy(wrf_ncfile, lat_top_right_SALLJ, lon_top_right_SALLJ) 

    ############################# read in some variables #############################
    #print('reading in variables')
    # read in the variables subsetted to a specified area
    pres_subset = getvar(wrf_ncfile, 'pressure')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]
    hght_subset = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in m NOTE: this is AGL not MSL!!
    radar_refl_max = getvar(wrf_ncfile, 'REFD_MAX')[int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]

    speed, drct = getvar(wrf_ncfile, 'wspd_wdir', units='kt') # in kts

    speed_subset = speed[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])] # in kts

    drct_subset = drct[:,int(bottom_left_xy[1]):int(top_right_xy[1]),int(bottom_left_xy[0]):int(top_right_xy[0])]

#                # get wind speed and direction at chosen height 1 and calculate wind shear
#                wspd_chosen1 = interplevel(speed_subset, hght_subset, chosen_height1)
#                wdir_chosen1 = interplevel(drct_subset, hght_subset, chosen_height1)
#                
#                udiff = wspd_chosen1*np.cos(np.radians(270-wdir_chosen1)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
#                vdiff = wspd_chosen1*np.sin(np.radians(270-wdir_chosen1)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
#
#                shear_chosen1 = np.sqrt(udiff**2 + vdiff**2)
#                
#                print('wspd_chosen1', wspd_chosen1[30, 30].values)
#                print('wdir_chosen1', wdir_chosen1[30, 30].values)
#                print('wspd_sfc', speed_subset[3, 30, 30].values)
#                print('wdir_sfc', drct_subset[3, 30, 30].values)
#                print('shear_chosen1', shear_chosen1[30, 30].values)
#                
#                # get wind speed and direction at chosen height 2 and calculate wind shear
#                wspd_chosen2 = interplevel(speed_subset, hght_subset, chosen_height2)
#                wdir_chosen2 = interplevel(drct_subset, hght_subset, chosen_height2)
#                
#                udiff = wspd_chosen2*np.cos(np.radians(270-wdir_chosen2)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
#                vdiff = wspd_chosen2*np.sin(np.radians(270-wdir_chosen2)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
#
#                shear_chosen2 = np.sqrt(udiff**2 + vdiff**2)
#                
#                # get wind speed and direction at chosen height 3 and calculate wind shear
#                wspd_chosen3 = interplevel(speed_subset, hght_subset, chosen_height3)
#                wdir_chosen3 = interplevel(drct_subset, hght_subset, chosen_height3)
#                
#                udiff = wspd_chosen3*np.cos(np.radians(270-wdir_chosen3)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
#                vdiff = wspd_chosen3*np.sin(np.radians(270-wdir_chosen3)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
#
#                shear_chosen3 = np.sqrt(udiff**2 + vdiff**2)
#                
#                # get wind speed and direction at chosen height 6 and calculate wind shear
#                wspd_chosen6 = interplevel(speed_subset, hght_subset, chosen_height6)
#                wdir_chosen6 = interplevel(drct_subset, hght_subset, chosen_height6)
#                
#                udiff = wspd_chosen6*np.cos(np.radians(270-wdir_chosen6)) - speed_subset[3, :, :]*np.cos(np.radians(270-drct_subset[3, :, :]))
#                vdiff = wspd_chosen6*np.sin(np.radians(270-wdir_chosen6)) - speed_subset[3, :, :]*np.sin(np.radians(270-drct_subset[3, :, :]))
#
#                shear_chosen6 = np.sqrt(udiff**2 + vdiff**2)
#                
#                # get wind shear 2-6 km
#                
#                udiff = wspd_chosen6*np.cos(np.radians(270-wdir_chosen6)) - wspd_chosen2*np.cos(np.radians(270-wdir_chosen2))
#                vdiff = wspd_chosen6*np.sin(np.radians(270-wdir_chosen6)) - wspd_chosen2*np.sin(np.radians(270-wdir_chosen2))
#
#                shear_chosen2_6 = np.sqrt(udiff**2 + vdiff**2)

    ############################# find SALLJ #############################

    # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind minimum 
    interp_levels_min = np.arange(0,min_search_hgt+250,250)

    pres_below_min = interplevel(pres_subset, hght_subset, interp_levels_min)
    hght_below_min = interplevel(hght_subset, hght_subset, interp_levels_min)
    speed_below_min = interplevel(speed_subset, hght_subset, interp_levels_min)
    drct_below_min = interplevel(drct_subset, hght_subset, interp_levels_min)

    # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind maximum
    interp_levels_max = np.arange(0,max_search_hgt+250,250)

    pres_below_max = interplevel(pres_subset, hght_subset, interp_levels_max)
    hght_below_max = interplevel(hght_subset, hght_subset, interp_levels_max)
    speed_below_max = interplevel(speed_subset, hght_subset, interp_levels_max)
    drct_below_max = interplevel(drct_subset, hght_subset, interp_levels_max)

    #print(pres_below_max.shape)
    #print(hght_below_max.shape)

    #print('pres_below_max', pres_below_max.values)
    #print('hght_below_max', hght_below_max.values)
    #print('speed_below_max', speed_below_max.values)

    ################ with xarray ###################

    # get max wind
    max_wind = speed_below_max.max(dim='level')

    # get height of max wind
    level_max_wind = hght_below_max.isel(level=speed_below_max.argmax('level'))
#                level_max_wind = speed_below_max.idxmax('level')
#                level_max_wind = speed_below_max[np.where(np.equal(speed_below_max,max_wind))].get_index('level')
    #level_max_wind2 = speed_below_max.idxmax(dim='level')

    # get pressure at max wind
    pres_max_wind = pres_below_max.isel(level=speed_below_max.argmax('level'))

    # get direction at max wind
    drct_at_max_wind = drct_below_max.isel(level=speed_below_max.argmax('level'))

    #print('drct_at_max_wind', drct_at_max_wind)

    # set dim of level of max wind array to 'level' (height)
    level_max_wind_subtract = xarray.DataArray(level_max_wind, dims=['south_north', 'west_east'])

    #print('level_max_wind_subtract', level_max_wind_subtract)

    #print('hght_below_max', hght_below_max.values)

    # subtract the heights of max wind from the heights
    hght_minus_level_max_wind = hght_below_min - level_max_wind_subtract

    #print('hght_minus_level_max_wind', hght_minus_level_max_wind)

    speed_below_min_masked_below_max = speed_below_min.where(hght_minus_level_max_wind > 0., np.nan)
    hght_below_min_masked_below_max = hght_below_min.where(hght_minus_level_max_wind > 0., np.nan)

    speed_below_min_masked_above_max = speed_below_min.where(hght_minus_level_max_wind < 0., np.nan)
    hght_below_min_masked_above_max = hght_below_min.where(hght_minus_level_max_wind < 0., np.nan)

    #print('speed_below_min_masked_below_max', speed_below_min_masked_below_max)

    # get min wind above max
    min_wind = speed_below_min_masked_below_max.min(dim='level')

    min_wind_below_max = speed_below_min_masked_above_max.min(dim='level')

    # get index of min above max

    #min_wind_index = speed_below_min_masked_below_max.idxmin(dim='level')

    #min_wind_index_below_max = speed_below_min_masked_above_max.idxmin(dim='level')

    level_min_wind = hght_below_min_masked_below_max.isel(level=speed_below_min_masked_below_max.argmin('level'))

    #level_min_wind_below_max = hght_below_min_masked_above_max.isel(level=speed_below_min_masked_above_max.argmin('level'))

    #print('min_wind', min_wind.values)

    # checks if max_wind meets threshold and keeps value if it does meet the threshold
    max_wind_meeting_threshold = max_wind.where(max_wind > max_wind_threshold, np.nan)

    #print('max_wind_meeting_threshold', max_wind_meeting_threshold.values)

    #print('max_wind_meeting_threshold', max_wind_meeting_threshold)

    # calculates decrease to min wind
    decrease_to_min = max_wind_meeting_threshold - min_wind

    decrease_to_min_below_max = max_wind_meeting_threshold - min_wind_below_max

    #print('decrease_to_min', decrease_to_min)

    # checks if decrease_to_min meets threshold and keeps value if it does meet the threshold
    decrease_to_min_meeting_threshold = decrease_to_min.where(decrease_to_min > decrease_to_min_threshold, np.nan)

    #print('decrease_to_min_meeting_threshold', decrease_to_min_meeting_threshold)

    # checks to see if the values met the other criteria (was not replaced with nan) and if it does leave the value
    drct_at_max_wind_meeting_threshold = drct_at_max_wind.where(np.isnan(decrease_to_min_meeting_threshold) == False, np.nan)

    #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold.values)

    # checks to see if wind at max_wind is from a northerly direction and keeps value if it is
    drct_at_max_wind_meeting_threshold = drct_at_max_wind_meeting_threshold.where((drct_at_max_wind_meeting_threshold <= 45) | (drct_at_max_wind_meeting_threshold >= 315), np.nan)

    #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold)
    
    num_points_SALLJ = np.count_nonzero(~np.isnan(drct_at_max_wind_meeting_threshold))
    num_total_points = drct_at_max_wind_meeting_threshold.size
    
    #print('num points', num_total_points)
    
    proporation_SALLJ = num_points_SALLJ/num_total_points
    
    #print('proporation coverage: ', proporation_SALLJ)
    
    prop_area_SALLJ[count] = proporation_SALLJ
    
    if proporation_SALLJ >= 0.12:
        
        print('found SALLJ')
        
        SALLJ_count = SALLJ_count + 1
        
        ### get SALLJ characteristics for points where the SALLJ criteria is met

        # get max_wind for points that meet SALLJ criteria
        max_wind_SALLJs = max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
        
        median_SALLJ_max_wind[count] = np.nanmedian(max_wind_SALLJs)
        
        # get level of max_wind for points that meet SALLJ criteria
        level_max_wind_SALLJs = level_max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
        
        median_SALLJ_height[count] = np.nanmedian(level_max_wind_SALLJs)
        
        #print('Median SALLJ wind: %s   Median SALLJ level: %s' %(median_SALLJ_max_wind, median_SALLJ_level))
        
        if proporation_SALLJ < 0.3:
            
            low_cov_SALLJ_count = low_cov_SALLJ_count + 1
            
            print('low_coverage_SALLJ')
        elif proporation_SALLJ <= 0.3 and proporation_SALLJ < 0.6:
            
            print('med_coverage_SALLJ')
            
            med_cov_SALLJ_count = med_cov_SALLJ_count + 1
            
        elif proporation_SALLJ >= 0.6:
            
            print('high_coverage_SALLJ')
            
            high_cov_SALLJ_count = high_cov_SALLJ_count + 1
        
    else:
        print('no SALLJ found')
        
#        SALLJ_status = 'no_SALLJ'
        #pass
        
        
        
    # calculate shear characteristics
    
    # calculate Richardson number
    
    # save MCS characteristics (time, centroid, growth rate, cloud sheild area) and environmental characteristics (SALLJ characterisitcs, shear) to netCDF or dateframe file
    
    wrf_ncfile.close()
    
print('SALLJ_count', SALLJ_count)
print('low_cov_SALLJ_count', low_cov_SALLJ_count)
print('med_cov_SALLJ_count', med_cov_SALLJ_count)
print('high_cov_SALLJ_count', high_cov_SALLJ_count)
    
pickle.dump(MCS_duration_filtered, open("MCS_duration2%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_majoraxislength_growth_filtered, open("MCS_majoraxislength_growth%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(mean_bulk_shear_0_3km, open("MCS_mean_bulk_shear_0_3km%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(mean_bulk_shear_0_6km, open("MCS_mean_bulk_shear_0_6km%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(mean_bulk_shear_2_6km, open("MCS_mean_bulk_shear_2_6km%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(prop_area_SALLJ, open("MCS_prop_area_SALLJ%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(median_SALLJ_max_wind, open("MCS_median_SALLJ_max_wind%s_%s.dat" %(filter_label, MCS_init_area), "wb"))
pickle.dump(median_SALLJ_height, open("MCS_median_SALLJ_height%s_%s.dat" %(filter_label, MCS_init_area), "wb"))

    