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

# NOTE: some variables may not be available script that created them was before they were added. Can either rerun that script or comment out certain varibale 

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
MCS_majoraxislength_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_0_3km = pickle.load(open(general_path + specific_inpath + "%s_mean_bulk_shear_0_3km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_0_6km = pickle.load(open(general_path + specific_inpath + "%s_mean_bulk_shear_0_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
mean_bulk_shear_2_6km = pickle.load( open(general_path + specific_inpath + "%s_mean_bulk_shear_2_6km%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_height = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850 = pickle.load( open(general_path + specific_inpath + "%s_q_850%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

fig, ax = plt.subplots()

#------ MCS Length Growth vs 850hPa q -------

#x_var = MCS_q_850
#y_var = MCS_majoraxislength_growth_filtered
#
#filename_x = 'MCS_q_850'
#filename_y = 'MCS_majoraxislength_growth_filtered'
#
#xlabel = 'mean q at 850hPa'
#ylabel = 'MCS_majoraxislength_growth (km/h)'

#------ MCS Area Growth vs 850hPa q -------

x_var = MCS_q_850
y_var = MCS_ccs_area_growth_filtered

filename_x = 'MCS_q_850'
filename_y = 'MCS_ccs_area_growth_filtered'

xlabel = 'mean q at 850hPa'
ylabel = 'MCS_ccs_area_growth_filtered (km2/h)'

#------ MCS Area Growth 2hr vs 850hPa q -------

#x_var = MCS_q_850
#y_var = MCS_ccs_area_growth_filtered_2hr
#
#filename_x = 'MCS_q_850'
#filename_y = 'MCS_ccs_area_growth_filtered_2hr'
#
#xlabel = 'mean q at 850hPa'
#ylabel = 'MCS_ccs_area_growth_filtered_2hr (km2/2h)'

#------ MCS Length Growth vs SALLJ Coverage -------

#x_var = MCS_prop_area_SALLJ
#y_var = MCS_majoraxislength_growth_filtered
#
#filename_x = 'MCS_prop_area_SALLJ'
#filename_y = 'MCS_majoraxislength_growth_filtered'
#
#xlabel = 'Coverage proportion w/ SALLJ'
#ylabel = 'MCS_majoraxislength_growth (km/h)'

#------ MCS Area Growth 2hr vs SALLJ Coverage -------

#x_var = MCS_prop_area_SALLJ
#y_var = MCS_ccs_area_growth_filtered_2hr
#
#filename_x = 'MCS_prop_area_SALLJ'
#filename_y = 'MCS_ccs_area_growth_filtered_2hr'
#
#xlabel = 'Coverage proportion w/ SALLJ'
#ylabel = 'MCS_ccs_area_growth_filtered_2hr (km2/2h)'

#------SALLJ Coverage vs Shear Strength -------

#x_var = mean_bulk_shear_2_6km
#y_var = MCS_prop_area_SALLJ
#
#filename_x = 'mean_bulk_shear_2_6km'
#filename_y = 'MCS_prop_area_SALLJ'
#
#xlabel = 'Mean 2-6 km shear (kts)'
#ylabel = 'Coverage proportion w/ SALLJ'

#------ MCS Length Growth vs Shear Strength -------

#x_var = mean_bulk_shear_2_6km
#y_var = MCS_majoraxislength_growth_filtered
#
#filename_x = 'mean_bulk_shear_2_6km'
#filename_y = 'MCS_majoraxislength_growth_filtered'
#
#xlabel = 'Mean 2-6 km shear (kts)'
#ylabel = 'MCS_majoraxislength_growth (km/h)'

#------ MCS Duration vs Shear Strength -------

#x_var = mean_bulk_shear_0_3km
#y_var = MCS_duration_filtered
#
#filename_x = 'mean_bulk_shear_0_3km'
#filename_y = 'MCS_duration_filtered'
#
#xlabel = 'Mean 0-3 km shear (kts)'
#ylabel = 'MCS duration (h)'

#------ MCS Area Growth 2hr vs Shear Strength -------

#x_var = mean_bulk_shear_2_6km
#y_var = MCS_ccs_area_growth_filtered_2hr
#
#filename_x = 'mean_bulk_shear_2_6km'
#filename_y = 'MCS_ccs_area_growth_filtered_2hr'
#
#xlabel = 'Mean 2-6 km shear (kts)'
#ylabel = 'MCS_ccs_area_growth_filtered_2hr (km2/2h)'

#------ Shear Strength vs SALLJ Stength (for only SALLJ times) ------- # NOTE: must turn off second scatter where c=median_SALLJ_max_wind below and change first scatter to c=MCS_prop_area_SALLJ_onlySALLJ

#indices_SALLJ = np.where(~np.isnan(median_SALLJ_max_wind))
#
#mean_bulk_shear_2_6km_onlySALLJ = mean_bulk_shear_2_6km[indices_SALLJ]
#median_SALLJ_max_wind_onlySALLJ = median_SALLJ_max_wind[indices_SALLJ]
#MCS_prop_area_SALLJ_onlySALLJ = MCS_prop_area_SALLJ[indices_SALLJ]
#
#x_var = median_SALLJ_max_wind_onlySALLJ
#y_var = mean_bulk_shear_2_6km_onlySALLJ
#
#filename_x = 'median_SALLJ_max_wind_onlySALLJ'
#filename_y = 'mean_bulk_shear_2_6km_onlySALLJ'
#
#xlabel = 'median_SALLJ_max_wind (kts)'
#ylabel = 'Mean 2-6 km shear (kts)'
#
#ax.set_aspect(0.5)

############################## Color by SALLJ coverage ###############################

plotted_fig = ax.scatter(x_var,y_var, c=MCS_prop_area_SALLJ, cmap='Reds', zorder=2) # MCS_prop_area_SALLJ_onlySALLJ

cbar = fig.colorbar(plotted_fig)
cbar.ax.set_ylabel('Proportion SALLJ')


############################## Outline color by SALLJ strength ###############################
#cc = []
#
#for SALLJ_max in median_SALLJ_max_wind:
#    if SALLJ_max >= 20 and SALLJ_max < 30:
#        cc.append('green')
#    elif SALLJ_max >= 30 and SALLJ_max < 40:
#        cc.append('orange')
#    elif SALLJ_max >= 40:
#        cc.append('red')
#    else:
#        cc.append('black')

ax.scatter(x_var, y_var, marker='o', s=90, c=median_SALLJ_max_wind, cmap='Blues', zorder=1)

#######################  polynomial fit ###################### 
## Fit a polynomial of degree 3
#degree = 1
#coeffs = np.polyfit(x_var, y_var, degree)
#
## Create a polynomial function from the coefficients
#poly_fit = np.poly1d(coeffs)
#
## Generate values for the fitted curve
#x_fit = np.linspace(min(x_var), max(x_var), 100)
#y_fit = poly_fit(x_fit)

####################### linear regression ###################### 

# Perform linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x_var, y_var)

# R-squared value
r_squared = r_value**2

# Generate the regression line
x_regress = np.linspace(min(x_var), max(x_var), 100)
y_regress = slope * x_regress + intercept

# Plot fitted line
plt.plot(x_regress, y_regress, color='red', linewidth=2)

print('r_squared', r_squared)

plt.text(0.1, 0.9, f'$R^2 = {r_squared:.4f}$')

ax.set_ylabel(ylabel)
ax.set_xlabel(xlabel)

plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

plt.savefig(general_path + specific_outpath + '%s_vs_%s_scatter%s.png' %(filename_y, filename_x, filter_label), dpi=200)

print('saved')
    