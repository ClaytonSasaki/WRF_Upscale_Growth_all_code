#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Oct 8 2024

@author: crs326

Plots distributions of MCS environment variables for *all points* at MCS initation times seperated by rapid growth and slow growth MCS

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
env_search_area = '2.00fromMCScentroid' # '0.75fromMCScentroid', '2.00fromMCScentroid'

hours_offset = -3
offset_label = '_%dhrPrior' %(abs(hours_offset)) # '' or '_%dhrPrior' %(abs(hours_offset))

##############################################################################

# get corresponding file labels usinf chosen inputs
env_search_text = '__EnvArea_%s' %(env_search_area)

if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')

# get input files paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/MCS_characteristics_vs_environment_MCS_centroid/'

specific_inpath = '%sarea_%s%s/data/' %(MCS_file_label, MCS_init_area, env_search_text)

# read in files
MCS_duration_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_duration2%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_majoraxislength_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_majoraxislength_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_majoraxislength_growth_filtered_hr = pickle.load(open(general_path + specific_inpath + "%s%s_majoraxislength_growth_hr%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_majoraxislength_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s%s_majoraxislength_growth_2hr%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_meantb_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_meantb_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_meantb_growth_filtered_hr = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_meantb_growth_hr%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_meantb_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_meantb_growth_2hr%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_area_growth%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_area_growth_filtered_hr = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_area_growth_hr%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MCS_ccs_area_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s%s_ccs_area_growth_2hr%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
MUCAPE = pickle.load(open(general_path + specific_inpath + "%s%s_MUCAPE_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
bulk_shear_0_3km = pickle.load(open(general_path + specific_inpath + "%s%s_bulk_shear_0_3km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
bulk_shear_0_6km = pickle.load(open(general_path + specific_inpath + "%s%s_bulk_shear_0_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
bulk_shear_2_6km = pickle.load(open(general_path + specific_inpath + "%s%s_bulk_shear_2_6km_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))
q_850 = pickle.load(open(general_path + specific_inpath + "%s%s_q_850_MCS_centroid%s_%s.dat" %(MCS_file_label, offset_label, filter_label, MCS_init_area), "rb"))

# replace nans in MUCAPE with 0
MUCAPE = np.nan_to_num(MUCAPE)

fig, ax = plt.subplots()

var_x = MCS_ccs_area_growth_filtered
var_x_name = 'MCS_ccs_area_growth_filtered'
var_y = MUCAPE
var_y_name = 'MUCAPE'

print(var_x_name, var_x)
print(var_y_name, var_y)

ax.hist2d(var_x, var_y, bins=20)

plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s/plots/' %(MCS_file_label, MCS_init_area, env_search_text)

plt.savefig(general_path + specific_outpath + '%s_vs_%s_2dhist%s%s.png' %(var_y_name, var_x_name, filter_label, offset_label), dpi=200)

print('saved')
    