#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on fri Sept 30 11:59:01 2022

@author: crs326

Compares SALLJ height (and speed) with 0-3km wind shear

"""

import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from scipy import stats

matplotlib.rcParams.update({'font.size': 18})

#################### VARIABLES TO CHANGE ##########################

filter_label = '_filtered_init_loc_and_start_type'

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'
MCS_init_area = 'large_area1'
SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'
env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

seperate_by_growth_rate = False

plot_rapid_growth_dist = False
plot_slow_growth_dist = False
plot_all_growth_dist = True

seperate_by_elevation = True

events_removed = True

##############################################################################

# get corresponding file labels usinf chosen inputs
SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
env_search_text = '__EnvArea_%s' %(env_search_area)

if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')

# get input files paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_MCS_characteristics_vs_environment/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# read in files
MCS_prop_area_SALLJ = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
MCS_duration_filtered = pickle.load(open(general_path + specific_inpath + "%s_duration2%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_ccs_area_filtered = pickle.load(open(general_path + specific_inpath + "2%s_ccs_area%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_majoraxislength_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_growth_stage_time_length_filtered = pickle.load(open(general_path + specific_inpath + "%s_growth_stage_time_length%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_ccs_area_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_0_3km = pickle.load(open(general_path + specific_inpath + "%s_mean_bulk_shear_0_3km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_0_6km = pickle.load(open(general_path + specific_inpath + "%s_mean_bulk_shear_0_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_2_6km = pickle.load( open(general_path + specific_inpath + "%s_mean_bulk_shear_2_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_height = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850 = pickle.load( open(general_path + specific_inpath + "%s_q_850%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

# get input files paths
centroid_elevation_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/firstStorms_InitialGrowthRate_shear_direction_MCS_centroids_plus0.5deg'

centroid_elevation_specific_inpath = '/data/'

offset_label = ''

# read in files
MCS_centroid_elevation = pickle.load(open(centroid_elevation_path + centroid_elevation_specific_inpath + "%s%s_MCS_centroid_elevation%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))

######## remove certain problematic events ########

if events_removed == True:

    events_bool = np.full(len(MCS_ccs_area_growth_filtered), True)

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
    
    events_bool = np.full(len(MCS_ccs_area_growth_filtered), True)
    
    events_removed_label = ''

MCS_centroid_elevation = np.array(MCS_centroid_elevation)[events_bool]
MCS_ccs_area_growth_filtered = np.array(MCS_ccs_area_growth_filtered)[events_bool]
MCS_growth_stage_time_length_filtered = np.array(MCS_growth_stage_time_length_filtered)[events_bool]

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
median_MCS_ccs_area_growth_rate_filtered = np.nanmedian(MCS_ccs_area_growth_rate_filtered)

if seperate_by_growth_rate == True:

    condition_MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_rate_filtered >= median_MCS_ccs_area_growth_rate_filtered

    mask1 = condition_MCS_ccs_area_growth_rate_filtered

    print('mask_rapid_growth', mask1)

    rapid_i = np.where(mask1 == True)

    print('rapid_i', rapid_i)

    #MCS_prop_area_SALLJ_mask1 = MCS_prop_area_SALLJ[mask1]
    #MCS_duration_filtered_mask1 = MCS_duration_filtered[mask1]
    #MCS_majoraxislength_growth_filtered_mask1 = MCS_majoraxislength_growth_filtered[mask1]
    #MCS_ccs_area_growth_filtered_mask1 = MCS_ccs_area_growth_filtered[mask1]
    MCS_ccs_area_growth_filtered_mask1 = MCS_ccs_area_growth_filtered[mask1]
    MCS_growth_stage_time_length_filtered_mask1 = MCS_growth_stage_time_length_filtered[mask1]
    MCS_ccs_area_growth_rate_filtered_mask1 = MCS_ccs_area_growth_rate_filtered[mask1]
    #bulk_shear_0_3km_mask1 = bulk_shear_0_3km[mask1]
    #bulk_shear_0_6km_mask1 = bulk_shear_0_6km[mask1]
    #bulk_shear_2_6km_mask1 = bulk_shear_2_6km[mask1]
    #median_SALLJ_max_wind_mask1 = median_SALLJ_max_wind[mask1]
    #median_SALLJ_height_mask1 = median_SALLJ_height[mask1]
    #MCS_q_850_mask1 = MCS_q_850[mask1]

    print('# rapid growth MCS', len(MCS_ccs_area_growth_rate_filtered_mask1))

    condition_MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_rate_filtered < median_MCS_ccs_area_growth_rate_filtered

    mask2 = condition_MCS_ccs_area_growth_rate_filtered

    #MCS_prop_area_SALLJ_mask2 = MCS_prop_area_SALLJ[mask2]
    #MCS_duration_filtered_mask2 = MCS_duration_filtered[mask2]
    #MCS_majoraxislength_growth_filtered_mask2 = MCS_majoraxislength_growth_filtered[mask2]
    #MCS_ccs_area_growth_filtered_mask2 = MCS_ccs_area_growth_filtered[mask2]
    MCS_ccs_area_growth_filtered_mask2 = MCS_ccs_area_growth_filtered[mask2]
    MCS_growth_stage_time_length_filtered_mask2 = MCS_growth_stage_time_length_filtered[mask2]
    MCS_ccs_area_growth_rate_filtered_mask2 = MCS_ccs_area_growth_rate_filtered[mask2]
    #bulk_shear_0_3km_mask2 = bulk_shear_0_3km[mask2]
    #bulk_shear_0_6km_mask2 = bulk_shear_0_6km[mask2]
    #bulk_shear_2_6km_mask2 = bulk_shear_2_6km[mask2]
    #median_SALLJ_max_wind_mask2 = median_SALLJ_max_wind[mask2]
    #median_SALLJ_height_mask2 = median_SALLJ_height[mask2]
    #MCS_q_850_mask2 = MCS_q_850[mask2]

    print('# slow growth MCS', len(MCS_ccs_area_growth_rate_filtered_mask2))

    data_all_growth = MCS_growth_stage_time_length_filtered
    data_rapid_growth = MCS_growth_stage_time_length_filtered_mask1
    data_slow_growth = MCS_growth_stage_time_length_filtered_mask2
    variable_name = 'MCS_growth_stage_time_length_filtered' # MCS_prop_area_SALLJ, MCS_0_3km_shear, MCS_q_850 
    x_label = 'growth time (hrs)' # proportion of area w/SALLJ, bulk 0-3 km shear, 850-hPa q

    fig, ax = plt.subplots()

    if plot_all_growth_dist == True:

        median_value = np.nanmedian(data_all_growth)
        dist_text_all_growth = '_all'   
        plt.hist(data_all_growth, bins=30, alpha=0.3, density=False, label='all MCS', color='gray')
        plt.axvline(median_value, color='gray', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

    else:
        dist_text_all_growth = ''

    if plot_rapid_growth_dist == True:

        median_value = np.nanmedian(data_rapid_growth)
        dist_text_rapid_growth = '_rapid'   
        plt.hist(data_rapid_growth, bins=30, alpha=0.3, density=False, label='rapid growth MCS', color='blue')
        plt.axvline(median_value, color='blue', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

    else:
        dist_text_rapid_growth = ''

    if plot_slow_growth_dist == True:

        median_value = np.nanmedian(data_slow_growth)
        dist_text_slow_growth = '_slow'
        plt.hist(data_slow_growth, bins=30, alpha=0.3, density=False, label='slow growth MCS', color='red')
        plt.axvline(median_value, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

    else:
        dist_text_slow_growth = ''
        
elif seperate_by_elevation == True:
    
    var = MCS_growth_stage_time_length_filtered
    variable_name = 'MCS_growth_stage_time_length_filtered'
    x_label = 'growth time (hrs)' # 'growth time (hrs)', 'MCS ccs area growth rate (km^2/hr)'
    
    fig, ax = plt.subplots()
    
    condition_elevation = MCS_centroid_elevation < 450

    mask = condition_elevation

    var_less500m = np.array(var)[mask]
    
    plt.hist(var_less500m, bins=30, alpha=0.3, density=False, label='<500 m', color='#0504aa')
    
    median_value = np.nanmedian(var_less500m)
    plt.axvline(median_value, color='#0504aa', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')
    
    condition_elevation = MCS_centroid_elevation >= 450

    mask = condition_elevation

    var_grt500m = np.array(var)[mask]
    
    plt.hist(var_grt500m, bins=30, alpha=0.3, density=False, label='>=500 m', color='red')
    
    median_value = np.nanmedian(var_grt500m)
    plt.axvline(median_value, color='red', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

# Adding labels and title
plt.xlabel(x_label)
#plt.ylabel('Probability Density')
plt.ylabel('Frequency')

# Show legend
plt.legend()

plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

if seperate_by_growth_rate == True:

    plt.savefig(general_path + specific_outpath + '5%s_hist%s%s%s_growth%s.png' %(variable_name, dist_text_all_growth, dist_text_rapid_growth, dist_text_slow_growth, filter_label), dpi=200)
    
elif seperate_by_elevation == True:
    
    plt.savefig(general_path + specific_outpath + '5%s_hist_all_growth_by_elevation%s.png' %(variable_name, filter_label), dpi=200)

print('saved')
    