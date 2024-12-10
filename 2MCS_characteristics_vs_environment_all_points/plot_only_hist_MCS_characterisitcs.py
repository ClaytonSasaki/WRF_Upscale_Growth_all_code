#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 14 2024

@author: crs326

Plots distributions of MCS characteristics for *all points* at MCS initation times.

"""

import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
from scipy import stats

matplotlib.rcParams.update({'font.size': 15})

#################### VARIABLES TO CHANGE ##########################

filter_label = '_filtered_init_loc_and_start_type'

MCS_type = 'all_MCS' # 'all_MCS' or 'robustMCS'
MCS_init_area = 'large_area1'
SALLJ_search_area = '2deg4degOffset1degNFromCentroid' # '2deg4degOffset1degNFromCentroid', '1deg3degBottomCentroid', '60-65W28-30SFixed'
env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

plot_rapid_growth_dist = True
plot_slow_growth_dist = True
plot_all_growth_dist = False

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
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCS_characteristics_vs_environment_all_points/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# read in files
MCS_prop_area_SALLJ = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
MCS_duration_filtered = pickle.load(open(general_path + specific_inpath + "%s_duration2_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_majoraxislength_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_3km = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_3km_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_6km = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_6km_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_2_6km = pickle.load( open(general_path + specific_inpath + "%s_bulk_shear_2_6km_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_height = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850 = pickle.load( open(general_path + specific_inpath + "%s_q_850_all_points%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

variable_name = 'ccs_area_growth' # MCS_prop_area_SALLJ, MCS_0_3km_shear, MCS_q_850 
x_label = 'ccs area growth' # proportion of area w/SALLJ, bulk 0-3 km shear, 850-hPa q

fig, ax = plt.subplots()

median_value = np.nanmedian(MCS_ccs_area_growth_filtered) 
plt.hist(MCS_ccs_area_growth_filtered, bins=30, alpha=0.3, density=True, label='ccs area growth 30min', color='gray')
plt.axvline(median_value, color='gray', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

median_value = np.nanmedian(MCS_ccs_area_growth_filtered_2hr)   
plt.hist(MCS_ccs_area_growth_filtered_2hr, bins=30, alpha=0.3, density=True, label='ccs area growth 2hr', color='blue')
plt.axvline(median_value, color='blue', linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

# Adding labels and title
plt.xlabel(x_label)
plt.ylabel('Probability Density')

# Show legend
plt.legend()

plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

plt.savefig(general_path + specific_outpath + '%s_hist%s.png' %(variable_name, filter_label), dpi=200)

print('saved')
    