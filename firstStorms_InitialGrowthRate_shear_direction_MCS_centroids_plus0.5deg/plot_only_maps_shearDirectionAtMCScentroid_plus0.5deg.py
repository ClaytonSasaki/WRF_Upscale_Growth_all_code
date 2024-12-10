#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 24 2024

@author: crs326

Plots maps of shear direction at MCS centroid. Additionally seperates by MCS growth rate.

"""

import matplotlib
matplotlib.use('Agg') 

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pickle
from scipy import stats

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, interp1d, vinterp, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

matplotlib.rcParams.update({'font.size': 18})

#################### VARIABLES TO CHANGE ##########################

filter_label = '_filtered_init_loc_and_start_type'

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'
MCS_init_area = 'large_area1'

env_offset_MCS = False
hours_offset = -1

seperation_plot = False
growth_type = 'slow'

plot_dir_shear_label = False

shear_type = '0-3km'

events_removed = True

##############################################################################

if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')
    
if env_offset_MCS == True:
    offset_label = '_%dhrPrior' %(abs(hours_offset))
else:
    offset_label = ''
    
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
    

# get input files paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/3shear_direction_MCS_centroids_plus0.5deg'

specific_inpath = '/data/'

# read in files
MCS_center_lons_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_center_lats_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_area_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_growth_stage_time_length_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_growth_stage_time_length%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
dir_shear_0_3km = pickle.load(open(general_path + specific_inpath + "%s%s_shearDir_0_3km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
dir_shear_0_6km = pickle.load(open(general_path + specific_inpath + "%s%s_shearDir_0_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
dir_shear_2_6km = pickle.load(open(general_path + specific_inpath + "%s%s_shearDir_2_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))

######## remove certain problematic events ########

if events_removed == True:

    events_bool = np.full(len(dir_shear_0_3km), True)

    # 2
    #remove_events_nums = [72,73,81,84,88,93,95,99,136,141,142,153]

    # 3
    #remove_events_nums = [52,53,69,72,73,81,82,84,88,90,93,95,9799,101,104,105,106,107,127,131,136,138,141,142,143,148,153,159]

    # 4
    #remove_events_nums = [1,7,10,154,72,73,81,84,88,93,95,99,136,141,142,153]
    
    # 5
    remove_events_nums = [1,7,10,154,72,73,81,84,88,93,95,99,136,141,142,153,57]

    events_bool[remove_events_nums] = False
    
    events_removed_label = '_events_removed'
    
else:
    
    events_bool = np.full(len(dir_shear_0_3km), True)
    
    events_removed_label = ''

MCS_center_lons_initiation_filtered = np.array(MCS_center_lons_initiation_filtered)[events_bool]
MCS_center_lats_initiation_filtered = np.array(MCS_center_lats_initiation_filtered)[events_bool]
MCS_ccs_area_growth_filtered = np.array(MCS_ccs_area_growth_filtered)[events_bool]
MCS_growth_stage_time_length_filtered = np.array(MCS_growth_stage_time_length_filtered)[events_bool]
dir_shear_0_3km = np.array(dir_shear_0_3km)[events_bool]
dir_shear_0_6km = np.array(dir_shear_0_6km)[events_bool]
dir_shear_2_6km = np.array(dir_shear_2_6km)[events_bool]

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
median_MCS_ccs_area_growth_rate_filtered = np.nanmedian(MCS_ccs_area_growth_rate_filtered)

if shear_type == '0-3km':
    dir_shear_chosen_layer = dir_shear_0_3km
elif shear_type == '0-6km':
    dir_shear_chosen_layer = dir_shear_0_6km
elif shear_type == '2-6km':
    dir_shear_chosen_layer = dir_shear_2_6km
else:
    print('Please enter a valid layer')

############ Plot terrain ############################# 

# get the file path to plot terrain
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'
wrf_file_path = path_wrf + '20181101/wrfout_d01_2018-11-01_00:00:00' # date doesn't matter as just plotting terrain

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

axs.scatter(-64.212, -31.298, s=240, color='red', transform=crs.PlateCarree(), zorder=5)
axs.scatter(-63.726, -29.906, s=240, color='red', transform=crs.PlateCarree(), zorder=5)

wrf_ncfile.close()

### plot MCS centroids colored by wind shear direction ###

if seperation_plot == True:
    
    ############ Change variable for seperation ############

    var_seperation = MCS_ccs_area_growth_rate_filtered
    var_seperation_name = 'ccs_area_growth_rate'
    seperation_threshold = median_MCS_ccs_area_growth_rate_filtered

    ############### Seperate by variable ################
    
    var_seperation_label = '_%s_by_%s%d' %(growth_type, var_seperation_name, seperation_threshold)

    if growth_type == 'rapid':
        mask = var_seperation >= seperation_threshold
    elif growth_type == 'slow':
        mask = var_seperation < seperation_threshold
    else:
        print('Please enter valid growth type')

    MCS_center_lons_initiation_filtered_mask = MCS_center_lons_initiation_filtered[mask]
    MCS_center_lats_initiation_filtered_mask = MCS_center_lats_initiation_filtered[mask]
    dir_shear_chosen_layer_mask = dir_shear_chosen_layer[mask]
    
    ################ plot MCS init centroids colored by shear direction ####################

    cc_shear_dir = [
        'blue' if (shear_dir <= 45) or (shear_dir >= 300) else 'red' if (135 <= shear_dir <= 225) else 'black'
        for shear_dir in dir_shear_chosen_layer_mask
    ]
    axs.scatter(MCS_center_lons_initiation_filtered_mask, MCS_center_lats_initiation_filtered_mask, transform=crs.PlateCarree(), color=cc_shear_dir, marker='*', zorder=7, s=1600)
    
else: # seperation_plot == False
    
    MCS_center_lons_initiation_filtered_mask = MCS_center_lons_initiation_filtered[:]
    MCS_center_lats_initiation_filtered_mask = MCS_center_lats_initiation_filtered[:]
    dir_shear_chosen_layer_mask = dir_shear_chosen_layer[:]
    
    axs.scatter(MCS_center_lons_initiation_filtered_mask, MCS_center_lats_initiation_filtered_mask, transform=crs.PlateCarree(), color='black', marker='*', zorder=7, s=1600)
    
    var_seperation_label = ''

print('# of MCS', len(MCS_center_lons_initiation_filtered_mask))

############### plot shear direction for each MCS ####################

if plot_dir_shear_label == True:

    for lon, lat, dir_shear in zip(MCS_center_lons_initiation_filtered_mask, MCS_center_lats_initiation_filtered_mask, dir_shear_chosen_layer_mask):
        plt.text(lon+0.15, lat, round(dir_shear), transform=crs.PlateCarree())
        
    dir_shear_label = '_dir_shear_label'
        
else:
    
    dir_shear_label = ''

############### plot number of MCS ####################

plt.text(.01, .99, 'n = %d' %(len(MCS_center_lons_initiation_filtered_mask)), ha='left', va='top', transform=axs.transAxes, fontsize=30)

general_outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/3shear_direction_MCS_centroids_plus0.5deg'

specific_outpath = '/data/'

plt.savefig(general_outpath + '/10%s%s%s_initCentroids_%s_%s%s%s_300.png' %(MCS_file_label, offset_label, var_seperation_label, shear_type, MCS_init_area, dir_shear_label, events_removed_label), dpi=600)

print('saved')

plt.close()
    