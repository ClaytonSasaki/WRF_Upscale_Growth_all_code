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

plot_rapid_growth_dist = False
plot_slow_growth_dist = False
plot_all_growth_dist = True

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
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/2MCS_characteristics_vs_environment/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# read in files
MCS_prop_area_SALLJ = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
MCS_duration_filtered = pickle.load(open(general_path + specific_inpath + "%s_duration2%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_filtered = pickle.load(open(general_path + specific_inpath + "%s_ccs_area%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_majoraxislength_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_0_3km = pickle.load(open(general_path + specific_inpath + "%s_mean_bulk_shear_0_3km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_0_6km = pickle.load(open(general_path + specific_inpath + "%s_mean_bulk_shear_0_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_2_6km = pickle.load( open(general_path + specific_inpath + "%s_mean_bulk_shear_2_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_height = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850 = pickle.load( open(general_path + specific_inpath + "%s_q_850%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

condition_MCS_ccs_area_growth_filtered_2hr = MCS_ccs_area_growth_filtered_2hr >= 10000

mask1 = condition_MCS_ccs_area_growth_filtered_2hr

#MCS_prop_area_SALLJ_mask1 = MCS_prop_area_SALLJ[mask1]
MCS_ccs_area_filtered_mask1 = MCS_ccs_area_filtered[mask1]
#MCS_duration_filtered_mask1 = MCS_duration_filtered[mask1]
#MCS_majoraxislength_growth_filtered_mask1 = MCS_majoraxislength_growth_filtered[mask1]
#MCS_ccs_area_growth_filtered_mask1 = MCS_ccs_area_growth_filtered[mask1]
MCS_ccs_area_growth_filtered_2hr_mask1 = MCS_ccs_area_growth_filtered_2hr[mask1]
#bulk_shear_0_3km_mask1 = bulk_shear_0_3km[mask1]
#bulk_shear_0_6km_mask1 = bulk_shear_0_6km[mask1]
#bulk_shear_2_6km_mask1 = bulk_shear_2_6km[mask1]
#median_SALLJ_max_wind_mask1 = median_SALLJ_max_wind[mask1]
#median_SALLJ_height_mask1 = median_SALLJ_height[mask1]
#MCS_q_850_mask1 = MCS_q_850[mask1]

print('# rapid growth MCS', len(MCS_ccs_area_growth_filtered_2hr_mask1))

condition_MCS_ccs_area_growth_filtered_2hr = MCS_ccs_area_growth_filtered_2hr < 10000

mask2 = condition_MCS_ccs_area_growth_filtered_2hr

#MCS_prop_area_SALLJ_mask2 = MCS_prop_area_SALLJ[mask2]
MCS_ccs_area_filtered_mask2 = MCS_ccs_area_filtered[mask2]
#MCS_duration_filtered_mask2 = MCS_duration_filtered[mask2]
#MCS_majoraxislength_growth_filtered_mask2 = MCS_majoraxislength_growth_filtered[mask2]
#MCS_ccs_area_growth_filtered_mask2 = MCS_ccs_area_growth_filtered[mask2]
MCS_ccs_area_growth_filtered_2hr_mask2 = MCS_ccs_area_growth_filtered_2hr[mask2]
#bulk_shear_0_3km_mask2 = bulk_shear_0_3km[mask2]
#bulk_shear_0_6km_mask2 = bulk_shear_0_6km[mask2]
#bulk_shear_2_6km_mask2 = bulk_shear_2_6km[mask2]
#median_SALLJ_max_wind_mask2 = median_SALLJ_max_wind[mask2]
#median_SALLJ_height_mask2 = median_SALLJ_height[mask2]
#MCS_q_850_mask2 = MCS_q_850[mask2]

print('# slow growth MCS', len(MCS_ccs_area_growth_filtered_2hr_mask2))

data_all_growth = MCS_ccs_area_filtered
data_rapid_growth = MCS_ccs_area_filtered_mask1
data_slow_growth = MCS_ccs_area_filtered_mask2
variable_name = 'MCS_ccs_area_filtered' # MCS_prop_area_SALLJ, MCS_0_3km_shear, MCS_q_850 
x_label = 'MCS ccs area (km^2)' # proportion of area w/SALLJ, bulk 0-3 km shear, 850-hPa q

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

# Adding labels and title
plt.xlabel(x_label)
#plt.ylabel('Probability Density')
plt.ylabel('Frequency')

# Show legend
plt.legend()

plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

plt.savefig(general_path + specific_outpath + '%s_hist%s%s%s_growth%s.png' %(variable_name, dist_text_all_growth, dist_text_rapid_growth, dist_text_slow_growth, filter_label), dpi=200)

print('saved')
    