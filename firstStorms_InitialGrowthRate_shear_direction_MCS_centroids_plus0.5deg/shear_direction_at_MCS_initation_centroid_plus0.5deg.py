#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 22 2024

@author: crs326

Gets MCS lat/lon centroid and times at initiation for MCSs in a defined region.

Calculates environmental characteristics (shear direction for *just* the *MCS initiation location centroid*. Plots on map for single point at the MCS initation time

based on file 'relaxed_criteria_SALLJ_coverage_characteristics_and_shear_at_MCS_initation_MCS_centroid_efficient_v3.py' with plotting from 'MCS_centroids_map.py'

which was based on file with similar name in MCS_characteristics_vs_environment_all_points directory (relaxed_criteria_SALLJ_coverage_characteristics_and_shear_at_MCS_initation_MCS_centroid_efficient_v3.py)

originally based on 'relaxed_criteria_SALLJ_spatial_areas_coverage_efficient_v3.py', 'Find_all_SALLJ_times_relaxed_efficient_avg_SALLJ_level.py', and 'SALLJ_spatial_areas_coverage.py'

"""

import matplotlib
matplotlib.use('Agg') 

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
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
                 cartopy_ylim, latlon_coords, interplevel, interp1d, vinterp, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

def get_wind_direction(u, v):
    """Calculates the wind direction in degrees from u and v components.

    Args:
        u (float): East-West wind component (positive is eastward)
        v (float): North-South wind component (positive is northward)

    Returns:
        float: Wind direction in degrees (0 is North, 90 is East, 180 is South, 270 is West)
    """

    # calculates the angle in radians between the positive x-axis and the line connecting the origin to the point (u, v), then converts the angle from radians to degrees, and finally adjusts the angle to represent the direction from which the wind is blowing (meteorological convention)
    direction = np.degrees(np.arctan2(u, v)) + 180
    
    # ensures the direction is within the range of 0 to 360 degrees
    return direction % 360

#################### VARIABLES TO CHANGE ##########################

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'

MCS_init_location_filter = True

MCS_start_type_filter = True

offset_MCS_and_conditions = False

hours_offset = 1

MCS_init_area = 'large_area1' # 'large_area1', 'large_area2', 'Zhang_area', 'SDC_area1'

#SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'
#
#env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

######### get MCS initation times and centroid locations #########

if MCS_type == 'all_MCS':
    MCS_tracked_file = 'mcs_tracks_pf_20181015_20190430.nc'
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_tracked_file = 'robust_mcs_tracks_20181015_20190430.nc'
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')

# read in MCS tracks
MCS_tracks_file_path = '/home/disk/monsoon/relampago/analysis//mcs_tracks/wrf/stats/%s' %(MCS_tracked_file)
MCS_tracks_ncfile = Dataset(MCS_tracks_file_path,'r')

MCS_status = np.array(MCS_tracks_ncfile.variables['mcs_status'])

### find indicies of when MCS is first found (when MCS_status is first 1 or 2) ###
    
# Initialize a boolean array of the same shape, filled with False
first_MCS_bol_array = np.zeros_like(MCS_status, dtype=bool)

# Initialize a boolean array of the same shape, filled with False
time_0_bol_array = np.zeros_like(MCS_status, dtype=bool)

# Initalize array with the length of the number of tracks
MCS_growth_stage_time_length = np.full(MCS_status.shape[0], np.nan)

# Find the first occurrence of `1` in each row
for i in range(MCS_status.shape[0]):
    
    ## Check if `1` or `2` (MCS or squall line) is in the row ##
    
    #print('MCS_status[i]', MCS_status[i])
    
    contains = np.isin(MCS_status[i], [1, 2]) # create condition if 1 or 2 is the in the row (MCS track)
    
    #print('contains', contains)
    
    if contains.any(): # check if 1 or 2 is the in the row (MCS track)
        
        # If MCS/Squall line found, get first tracked time (by setting first index to True)
        time_0_bol_array[i, 0] = True
        
        # Get the index where MCS/Squall line and set that index to True
        first_one_index = np.where(contains)[0][0]
        #print('np.where(contains)', np.where(contains))
        #print('first_one_index', first_one_index)
        first_MCS_bol_array[i, first_one_index] = True
        #print('first_MCS_bol_array[i, :]', first_MCS_bol_array[i, :])
        
        # The index is uses to calculate the duration as each time step is 30 min
        MCS_growth_stage_time_length[i] = first_one_index/2
        
## Possible alternative way without loops
#mask = (MCS_status == 1) | (MCS_status == 2)
#
## Find the index of the first occurrence of 1 or 2 in each row
#first_occurrence_indices = np.argmax(mask, axis=1)
#
## Create an array of False with the same shape as arr
#bool_array = np.zeros_like(MCS_status, dtype=bool)
#
## Set the first occurrence of 1 or 2 in each row to True using advanced indexing
#bool_array[np.arange(MCS_status.shape[0]), first_occurrence_indices] = True
#
## Make sure True is only set where there was a 1 or 2 in the original array
#bool_array &= mask
#
#print(bool_array)

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

# save MCS initation times and centroid locations of those within defined box
if MCS_init_area == 'large_area1':
    lon_min = -66.0
    lon_max = -57.0
    lat_min = -36.0
    lat_max = -28.0

elif MCS_init_area == 'large_area2':
    lon_min = -66.0
    lon_max = -57.0
    lat_min = -36.0
    lat_max = -29.5
    
elif MCS_init_area == 'Zhang_area':
    lon_min = -70.0
    lon_max = -59.0
    lat_min = -36.0
    lat_max = -26.5
    
elif MCS_init_area == 'SDC_area1':
    lon_min = -66.25
    lon_max = -61.1
    lat_min = -34.0
    lat_max = -29.5

else:
    print('Please add the matching lats and lons to search!')
        
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

ccs_area_init = np.array(MCS_tracks_ncfile.variables['ccs_area'])[:,0]
MCS_ccs_area_init = np.array(MCS_tracks_ncfile.variables['ccs_area'])[first_MCS_bol_array]

print('ccs_area_init', ccs_area_init)
print('MCS_ccs_area_init', MCS_ccs_area_init)

#print('first track area', np.array(MCS_tracks_ncfile.variables['ccs_area'])[0,:])
#print('first track status', np.array(MCS_tracks_ncfile.variables['mcs_status'])[0,:])

MCS_ccs_area_growth = MCS_ccs_area_init - ccs_area_init

MCS_ccs_area_growth_filtered = MCS_ccs_area_growth[mask]

print('MCS_ccs_area_growth_filtered', MCS_ccs_area_growth_filtered)

MCS_growth_stage_time_length_filtered = MCS_growth_stage_time_length[mask]

print('MCS_growth_stage_time_length_filtered', MCS_growth_stage_time_length_filtered)

print('rate', MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered)
print('median', np.nanmedian(MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered))

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
MCS_tracks_ncfile.close

######### Get area average environmental data centered on each MCS track chosen ######### 

#bulk_shear_0_3km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#bulk_shear_0_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#bulk_shear_2_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#q_850 = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
#MUCAPE = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)

terrain_plotted = False

shearDir_0_3km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
shearDir_0_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
shearDir_2_6km = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)
MCS_dt = []
MCS_centroid_elevation = np.full(len(MCS_datetime_initiation_str_filtered), np.nan)

# go through list of times/centroids for each MCS to get corresponding environmental conditions
for count, (MCS_datetime, MCS_center_lon, MCS_center_lat) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered)):

    # get the file necessary
    path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  
                                   
    MCS_datetime_dt = datetime.strptime(MCS_datetime, '%Y-%m-%d_%H:%M:')
    
    MCS_dt.append(MCS_datetime_dt)
    
    if offset_MCS_and_conditions == True:
        
        conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)
        
        if hours_offset < 0:
            offset_label = '_%dhrPrior' %(abs(hours_offset))
        elif hours_offset > 0:
            offset_label = '_%dhrAfter' %(abs(hours_offset))
        else:
            print('Enter valid hours_offset')
        
    else: #offset_MCS_and_conditions == False
        
        conditions_datetime_dt = MCS_datetime_dt
        offset_label = ''
    
    MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')

    wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))

    wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name
    
    print('path of MCS within region: ', wrf_file_path)

    # get the netCDF
    wrf_ncfile = Dataset(wrf_file_path,'r')
    
    if terrain_plotted == False:
        
        # Get the terrain height
        terrain = getvar(wrf_ncfile, 'HGT')

        # Get the latitude and longitude points
        lats, lons = latlon_coords(terrain)

        # Get the cartopy mapping object
        cart_proj = get_cartopy(terrain)

        area = 'broad_zoom'
        lat_bottom_left = -38.0
        lon_bottom_left = -70.0

        lat_top_right = -26.0
        lon_top_right = -54.0

        barb_interval = 35

        fig, axs = plt.subplots(1, 1, figsize=(32, 32), subplot_kw={'projection': cart_proj})


        bounds = GeoBounds(CoordPair(lat=lat_bottom_left, lon=lon_bottom_left),
                               CoordPair(lat=lat_top_right, lon=lon_top_right))

        axs.set_xlim(cartopy_xlim(terrain, geobounds=bounds))
        axs.set_ylim(cartopy_ylim(terrain, geobounds=bounds))

        ########### plot features, terrain, and gridlines ###########

        axs.add_feature(OCEAN, zorder=2)
        axs.add_feature(LAKES, alpha=0.5, zorder=2)
        axs.add_feature(RIVERS, zorder=2)
        axs.add_feature(BORDERS, edgecolor='gray', zorder=2)

        # Make  filled contours for the terrain
        terrain_plot = axs.contourf(to_np(lons), to_np(lats), to_np(terrain), levels=[0, 500, 1000], colors=['papayawhip', 'navajowhite','burlywood'], extend='max', transform=crs.PlateCarree(), zorder=1)

        # Add a color bar
        #                fig.subplots_adjust(right=0.8)
        #                cbar_ax1 = fig.add_axes([0.15, 0.05, 0.25, 0.05])
        #                fig.colorbar(terrain_plot,cax=cbar_ax1, label='Elevation (m)', orientation='horizontal')

        #cbar = plt.colorbar(terrain_plot, ax=axs, shrink=.2, orientation='horizontal')
        #cbar.set_label('Elevation (m)')

        # Add the gridlines
        gl = axs.gridlines(crs=crs.PlateCarree(), linewidth=1, color='gray', alpha=0.5, draw_labels=True, linestyle='--', zorder=6)
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = True
        gl.bottom_labels = True


        # plot MCS initation box
        axs.plot([lon_min, lon_max], [lat_max, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_min, lon_min], [lat_min, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_min, lon_max], [lat_min, lat_min], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_max, lon_max], [lat_min, lat_max], 'k-', lw=6, transform=crs.PlateCarree(), zorder=6)


        #######################################################

        #plt.title("%s" % (area))

        axs.scatter(-64.212, -31.298, s=240, color='red', transform=crs.PlateCarree(), zorder=5)
        axs.scatter(-63.726, -29.906, s=240, color='red', transform=crs.PlateCarree(), zorder=5)
        
        terrain_plotted = True
    
    # get xy for centroid
    centroid_xy = ll_to_xy(wrf_ncfile, MCS_center_lat+0.5, MCS_center_lon)
    
    # Get the terrain height at MCS centroid
    terrain_centroid = getvar(wrf_ncfile, 'HGT')[int(centroid_xy[1]),int(centroid_xy[0])].values
    
    MCS_centroid_elevation[count] = terrain_centroid

    # read in the variables 
    hght_original = getvar(wrf_ncfile, 'height', msl=False, units='m')[:,int(centroid_xy[1]),int(centroid_xy[0])].values # in m NOTE: this is AGL not MSL!!
    u_original = getvar(wrf_ncfile, 'ua', units='kt')[:,int(centroid_xy[1]),int(centroid_xy[0])].values  # in kts
    v_original = getvar(wrf_ncfile, 'va', units='kt')[:,int(centroid_xy[1]),int(centroid_xy[0])].values # in kts
    
    #print(hght_original.shape)
    
    #hght_original = np.transpose(hght_original)
    #u_original = np.transpose(u_original)
    #v_original = np.transpose(v_original)
    
    #print(v_original)
    
    #print(u_original)
    #print(hght_original)
    
    u2km = interp1d(u_original, hght_original, np.asarray([2000.])).values
    v2km = interp1d(v_original, hght_original, np.asarray([2000.])).values

    u3km = interp1d(u_original, hght_original, np.asarray([3000.])).values
    v3km = interp1d(v_original, hght_original, np.asarray([3000.])).values
    
    u6km = interp1d(u_original, hght_original, np.asarray([6000.])).values
    v6km = interp1d(v_original, hght_original, np.asarray([6000.])).values
    
    u_sfc = u_original[3]
    v_sfc = v_original[3]
    
    #print('v_sfc', v_sfc)
    
    # calculate 0-3 km shear at all points
    udiff_0_3 = u3km - u_sfc
    vdiff_0_3 = v3km - v_sfc
    
    dir_shear_0_3km = get_wind_direction(udiff_0_3, vdiff_0_3)
    
    shearDir_0_3km[count] = dir_shear_0_3km[0]
    
    # calculate 0-6 km shear at all points
    udiff_0_6 = u6km - u_sfc
    vdiff_0_6 = v6km - v_sfc
    
    dir_shear_0_6km = get_wind_direction(udiff_0_6, vdiff_0_6)
    
    shearDir_0_6km[count] = dir_shear_0_6km[0]
    
    # calculate 2-6 km shear at all points
    udiff_2_6 = u6km - u2km
    vdiff_2_6 = v6km - v2km
    
    dir_shear_2_6km = get_wind_direction(udiff_2_6, vdiff_2_6)
    
    shearDir_2_6km[count] = dir_shear_2_6km[0]
    
    print('dir_shear_0_3km', dir_shear_0_3km)
    
    if(dir_shear_0_3km <= 45) or dir_shear_0_3km >= 315:
        cc = 'blue'
    elif 135 <= dir_shear_0_3km <= 225:
        cc = 'red'
    else:
        cc = 'black'
        
    # plot MCS centroid
    axs.scatter(MCS_center_lon, MCS_center_lat, transform=crs.PlateCarree(), color=cc, marker='*', zorder=7, s=1600)
    
#    shear3km = np.sqrt(udiff_0_3**2 +vdiff_0_3**2)
#
#    bulk_shear_0_3km[count] = shear3km[0]
#    
#    print('shear3km', shear3km[0])
#    
#    udiff_0_6 = u6km - u_sfc
#    vdiff_0_6 = v6km - v_sfc
#    
#    shear6km = np.sqrt(udiff_0_6**2 + vdiff_0_6**2)
#
#    bulk_shear_0_6km[count] = shear6km
#    
#    print('shear6km', shear6km[0])
#    
#    udiff_2_6 = u6km - u2km
#    vdiff_2_6 = v6km - v2km
#    
#    shear2_6km = np.sqrt(udiff_2_6**2 + vdiff_2_6**2)
#    
#    print('shear2_6km', shear2_6km[0])
#
#    bulk_shear_2_6km[count] = shear2_6km
    
#    ######### calculate q at 850 hPa #########
#    level = np.asarray([850.])
#    
#    # read in the variables for the subsetted area 
#    pres = getvar(wrf_ncfile, 'pressure')[:,int(centroid_xy[1]),int(centroid_xy[0])].values
#
#    temp_c_subset = getvar(wrf_ncfile, 'temp', units='degC')[:,int(centroid_xy[1]),int(centroid_xy[0])].values
#
#    RH_subset = getvar(wrf_ncfile, 'rh')[:,int(centroid_xy[1]),int(centroid_xy[0])].values # in kg/kg
#
#    # interpolate to the set level
#    temp_c_subset_level = np.squeeze(interp1d(temp_c_subset, pres, level)).values
#
#    RH_subset_level = np.squeeze(interp1d(RH_subset, pres, level)).values
#
#    e_sat_subset_level = SatVap(temp_c_subset_level)
#
#    e_subset_level = (RH_subset_level/100.) * e_sat_subset_level
#
#    w_subset_level = MixRatio(e_subset_level,level[0]*100.)
#    #w = (e*Rs_da) / (Rs_v*(pres*100. - e))
#    q_subset_level = (w_subset_level / (w_subset_level+1))*1000 #g/kg
#
#    q_850[count] = q_subset_level
#    
#    print('q_subset_level', q_subset_level)
#    
#    MCAPE_original = getvar(wrf_ncfile, 'cape_2d')[0][int(centroid_xy[1]),int(centroid_xy[0])].values
#    
#    MUCAPE[count] = MCAPE_original
#    
#    print('MCAPE_original', MCAPE_original)

    wrf_ncfile.close()
    
general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

specific_outpath = '/data/'

plt.savefig(general_outpath + '/%s%s_initCentroids_%s.png' %(MCS_file_label, offset_label, MCS_init_area), dpi=600)

print('saved')

plt.close()

pickle.dump(MCS_dt, open(general_outpath + specific_outpath + "%s%s_MCS_dt%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
pickle.dump(MCS_centroid_elevation, open(general_outpath + specific_outpath + "%s%s_MCS_centroid_elevation%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(MCS_center_lons_initiation_filtered, open(general_outpath + specific_outpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(MCS_center_lats_initiation_filtered, open(general_outpath + specific_outpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(MCS_ccs_area_growth_filtered, open(general_outpath + specific_outpath + "%s%s_ccs_area_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(MCS_growth_stage_time_length_filtered, open(general_outpath + specific_outpath + "%s%s_growth_stage_time_length%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(shearDir_0_3km, open(general_outpath + specific_outpath + "%s%s_shearDir_0_3km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(shearDir_0_6km, open(general_outpath + specific_outpath + "%s%s_shearDir_0_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))
#pickle.dump(shearDir_2_6km, open(general_outpath + specific_outpath + "%s%s_shearDir_2_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "wb"))