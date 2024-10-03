#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Oct 10 15:00:00 2022

@author: crs326

Makes maps that plots wind shear at a chosen height and shows points that meet the SALLJ criteria and the height at which the criteria is met from the PNNL WRF. The same as shear_and_SALLJ_strength_map.py but with one close zoom plotted and updated for readability.

Plots maps for all chosen MCS initation times

Also plots location where MCS centroid was searched (black box), MCS centroid (black star), and where SALLJ coverage is calculated (dashed black box)

Based upon REFD_MAX_and_shear_map.py and SALLJ_spatial_map_zoom_efficient_v3.py
"""

import matplotlib
matplotlib.use('Agg') 

import math
import os

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib.colors as mc
import metpy.calc as mpcalc
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from numpy import exp, where, ma, cos, sin, pi, amax, amin
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d
#import time
import xarray

# ------------------------- cmaps ------------------------------

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

import matplotlib.colors as colors

#start_time = time.time()
matplotlib.rcParams.update({'font.size': 47})

# used for contourf of terrain 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap = plt.get_cmap('terrain')
new_cmap = truncate_colormap(cmap, 0.21, 1)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

# used for SALLJ
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('Greys')
new_cmap2 = truncate_colormap(cmap2, 0.2, 1.0)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# used for contour of wind
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('Greens')
new_cmap2 = truncate_colormap(cmap2, 0.1, 1.0)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

cdict1 = {'charteruse', 'lime', 'mediumsegreen', 'forestgreen', 'darkgreen'}
new_green_cmap = mc.LinearSegmentedColormap('newGreens', cdict1)

#################### VARIABLES TO CHANGE ##########################

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'

MCS_init_location_filter = True

MCS_start_type_filter = True

offset_MCS_and_conditions = False
hours_offset = -1

MCS_init_area = 'largeArea2'
SALLJ_search_area = '2deg4degOffset1degNFromCentroid'
plot_SALLJ_search_area = True

env_search_area = '0.75fromMCScentroid' # '0.75fromMCScentroid'
plot_env_search_area = True

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
print('lons of MCS init \n', MCS_center_lons_initiation)
print('lats of MCS init \n', MCS_center_lats_initiation)

######### create filter MCS tracks by centroid location ######### 

# save MCS initation times and centroid locations of those within defined box
if MCS_init_area == 'largeArea1':
    lon_min = -66.0
    lon_max = -57.0
    lat_min = -36.0
    lat_max = -28.0

elif MCS_init_area == 'largeArea2':
    lon_min = -66.0
    lon_max = -57.0
    lat_min = -36.0
    lat_max = -29.5
    
elif MCS_init_area == 'ZhangArea':
    lon_min = -70.0
    lon_max = -59.0
    lat_min = -36.0
    lat_max = -26.5

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

#print('MCS_duration', MCS_duration)
#print('len(MCS_duration)', len(MCS_duration))

MCS_duration_filtered = MCS_duration[mask]

MCS_majoraxislength_init = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,0]

MCS_majoraxislength_init_2 = np.array(MCS_tracks_ncfile.variables['majoraxislength'])[:,1]

MCS_majoraxislength_growth = MCS_majoraxislength_init_2 - MCS_majoraxislength_init

MCS_majoraxislength_growth_filtered = MCS_majoraxislength_growth[mask]

# close the MCS tracks file
MCS_tracks_ncfile.close()

# ----------------- plotting terrain from first file -----------------

######################################

# get the file necessary
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'                  

MCS_datetime_dt = datetime.strptime(MCS_datetime_initiation_str_filtered[0], '%Y-%m-%d_%H:%M:')

if offset_MCS_and_conditions == True:

    conditions_datetime_dt = MCS_datetime_dt + timedelta(hours = hours_offset)

else: #offset_MCS_and_conditions == False

    conditions_datetime_dt = MCS_datetime_dt

MCS_datetime_str = conditions_datetime_dt.strftime('%Y%m%d')

wrf_file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(conditions_datetime_dt.strftime('%Y'), conditions_datetime_dt.strftime('%m'), conditions_datetime_dt.strftime('%d'), conditions_datetime_dt.strftime('%H'))

wrf_file_path = path_wrf + MCS_datetime_str + wrf_file_name

#print('path of MCS within region: ', wrf_file_path)

# get the netCDF
wrf_ncfile = Dataset(wrf_file_path,'r')

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

wrf_ncfile.close()

print('terrain plotted')

# go through list of times/centroids for each MCS to get corresponding environmental conditions
for count, (MCS_datetime, MCS_center_lon, MCS_center_lat) in enumerate(zip(MCS_datetime_initiation_str_filtered, MCS_center_lons_initiation_filtered, MCS_center_lats_initiation_filtered)):
    
        
    # plot MCS centroid
    axs.scatter(MCS_center_lon, MCS_center_lat, transform=crs.PlateCarree(), color='black', marker='*', zorder=7, s=1600)
    
    # plot SALLJ search area 
    if plot_SALLJ_search_area == True:
        
        SALLJ_search_text = '_wSALLJarea_%s' %(SALLJ_search_area)
    
        if SALLJ_search_area == '1deg3degBottomCentroid':

            lat_bottom_left_SALLJ = MCS_center_lat
            lon_bottom_left_SALLJ = MCS_center_lon - 1.50
            lat_top_right_SALLJ = MCS_center_lat + 1.00 
            lon_top_right_SALLJ = MCS_center_lon + 1.50

        elif SALLJ_search_area == '2deg4degOffset1degNFromCentroid':

            lat_bottom_left_SALLJ = MCS_center_lat + 1.00
            lon_bottom_left_SALLJ = MCS_center_lon - 2.00
            lat_top_right_SALLJ = MCS_center_lat + 3.00 
            lon_top_right_SALLJ = MCS_center_lon + 2.00
        else:
            print('Please enter valid SALLJ search area') # will throw error
        
        
        # plot MCS initation box
        axs.plot([lon_bottom_left_SALLJ, lon_top_right_SALLJ], [lat_top_right_SALLJ, lat_top_right_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_SALLJ, lon_bottom_left_SALLJ], [lat_bottom_left_SALLJ, lat_top_right_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_SALLJ, lon_top_right_SALLJ], [lat_bottom_left_SALLJ, lat_bottom_left_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_top_right_SALLJ, lon_top_right_SALLJ], [lat_bottom_left_SALLJ, lat_top_right_SALLJ], 'k--', lw=6, transform=crs.PlateCarree(), zorder=6)
            
    else:
        
        SALLJ_search_text = ''
        
        
        
    # plot environmental search area    
    if plot_env_search_area == True:
        
        env_search_text = '_wEnvironmentArea_%s' %(env_search_area)
    
        if env_search_area == '0.75fromMCScentroid':

            # get lats/lons of region based on centroid
            lat_bottom_left_env = MCS_center_lat - 0.75 
            lon_bottom_left_env = MCS_center_lon - 0.75

            lat_top_right_env = MCS_center_lat + 0.75 
            lon_top_right_env = MCS_center_lon + 0.75

        else:
            print('Please enter valid env search area') # will throw error
        
        
        # plot MCS initation box
        axs.plot([lon_bottom_left_env, lon_top_right_env], [lat_top_right_env, lat_top_right_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_env, lon_bottom_left_env], [lat_bottom_left_env, lat_top_right_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_bottom_left_env, lon_top_right_env], [lat_bottom_left_env, lat_bottom_left_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
        axs.plot([lon_top_right_env, lon_top_right_env], [lat_bottom_left_env, lat_top_right_env], 'c--', lw=6, transform=crs.PlateCarree(), zorder=6)
            
    else:
        
        env_search_text = ''

    #print 'saving'

    #plt.suptitle("Max Reflectivity and sfc-%d wind shear: %sZ" % (level, file_name[11:24]))
    
outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCS_centroids_map'

plt.savefig(outpath + '/2%s_initCentroids_%s%s%s.png' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text), dpi=600)

print('saved')

plt.close()

