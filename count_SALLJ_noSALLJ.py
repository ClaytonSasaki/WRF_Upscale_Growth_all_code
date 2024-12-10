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

seperation_plot = True
growth_type = 'rapid'

plot_label = False

plot_var = 'SALLJ'

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
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/3shear_direction_MCS_centroids'

specific_inpath = '/data/'

# read in files
MCS_center_lons_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lons%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_center_lats_initiation_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_MCS_center_lats%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_area_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_growth_stage_time_length_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_growth_stage_time_length%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))

general_path_2 = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/3MCS_characteristics_vs_environment_all_points/'

SALLJ_search_area = '2deg4degOffset1degNFromCentroid'
env_search_area = '2.00fromMCScentroid'

# get corresponding file labels usinf chosen inputs
SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
env_search_text = '__EnvArea_%s' %(env_search_area)

specific_inpath_2 = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

MCS_prop_area_SALLJ_all_points_byEvent = pickle.load(open(general_path_2 + specific_inpath_2 + "%s_prop_area_SALLJ_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

# get first value for each event as they are all the same repeated
MCS_prop_area_SALLJ = np.array([arr[0,0] for arr in MCS_prop_area_SALLJ_all_points_byEvent])

######## remove certain problematic events ########

if events_removed == True:

    events_bool = np.full(len(MCS_prop_area_SALLJ), True)

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
    
    events_bool = np.full(len(MCS_prop_area_SALLJ), True)
    
    events_removed_label = ''

MCS_center_lons_initiation_filtered = np.array(MCS_center_lons_initiation_filtered)[events_bool]
MCS_center_lats_initiation_filtered = np.array(MCS_center_lats_initiation_filtered)[events_bool]
MCS_ccs_area_growth_filtered = np.array(MCS_ccs_area_growth_filtered)[events_bool]
MCS_growth_stage_time_length_filtered = np.array(MCS_growth_stage_time_length_filtered)[events_bool]
MCS_prop_area_SALLJ = np.array(MCS_prop_area_SALLJ)[events_bool]

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
median_MCS_ccs_area_growth_rate_filtered = np.nanmedian(MCS_ccs_area_growth_rate_filtered)

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
    MCS_prop_area_SALLJ_mask = MCS_prop_area_SALLJ[mask]
    
else: # seperation_plot == False
    
    MCS_center_lons_initiation_filtered_mask = MCS_center_lons_initiation_filtered[:]
    MCS_center_lats_initiation_filtered_mask = MCS_center_lats_initiation_filtered[:]
    MCS_prop_area_SALLJ_mask = MCS_prop_area_SALLJ[:]
    
    var_seperation_label = ''

print('# of MCS', len(MCS_center_lons_initiation_filtered_mask))


print((MCS_prop_area_SALLJ_mask > 0.2).sum())
    