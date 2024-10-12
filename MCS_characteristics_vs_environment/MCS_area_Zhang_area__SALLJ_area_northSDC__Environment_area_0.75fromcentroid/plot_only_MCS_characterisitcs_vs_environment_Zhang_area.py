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

filter_label = '_filtered_init_loc_and_start_type'

MCS_init_area = 'Zhang_area'

MCS_prop_area_SALLJ = pickle.load(open("MCS_prop_area_SALLJ%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
#MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
MCS_duration_filtered = pickle.load(open("MCS_duration2%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
MCS_majoraxislength_growth_filtered = pickle.load(open("MCS_majoraxislength_growth%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
mean_bulk_shear_0_3km = pickle.load(open("MCS_mean_bulk_shear_0_3km%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
mean_bulk_shear_0_6km = pickle.load(open("MCS_mean_bulk_shear_0_6km%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
mean_bulk_shear_2_6km = pickle.load( open("MCS_mean_bulk_shear_2_6km%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
median_SALLJ_max_wind = pickle.load( open("MCS_median_SALLJ_max_wind%s_%s.dat" %(filter_label, MCS_init_area), "rb"))
median_SALLJ_height = pickle.load( open("MCS_median_SALLJ_height%s_%s.dat" %(filter_label, MCS_init_area), "rb"))

fig, ax = plt.subplots()

#------ MCS Length Growth vs SALLJ Coverage -------

#x_var = MCS_prop_area_SALLJ
#y_var = MCS_majoraxislength_growth_filtered
#
#filename_x = 'MCS_prop_area_SALLJ'
#filename_y = 'MCS_majoraxislength_growth_filtered'
#
#xlabel = 'Coverage proportion w/ SALLJ'
#ylabel = 'MCS_majoraxislength_growth (km/h)'

#------SALLJ Coverage vs Shear Strength -------

#x_var = mean_bulk_shear_0_3km
#y_var = MCS_prop_area_SALLJ
#
#filename_x = 'mean_bulk_shear_0_3km'
#filename_y = 'MCS_prop_area_SALLJ'
#
#xlabel = 'Mean 0-3 km shear (kts)'
#ylabel = 'Coverage proportion w/ SALLJ'

#------ MCS Length Growth vs Shear Strength -------

#x_var = mean_bulk_shear_0_3km
#y_var = MCS_majoraxislength_growth_filtered
#
#filename_x = 'mean_bulk_shear_0_3km'
#filename_y = 'MCS_majoraxislength_growth_filtered'
#
#xlabel = 'Mean 0-3 km shear (kts)'
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

#------ Shear Strength vs SALLJ Stength (for only SALLJ times) ------- # NOTE: must turn off second scatter by median_SALLJ_max_wind below

indices_SALLJ = np.where(~np.isnan(median_SALLJ_max_wind))

mean_bulk_shear_2_6km_onlySALLJ = mean_bulk_shear_2_6km[indices_SALLJ]
median_SALLJ_max_wind_onlySALLJ = median_SALLJ_max_wind[indices_SALLJ]
MCS_prop_area_SALLJ_onlySALLJ = MCS_prop_area_SALLJ[indices_SALLJ]

x_var = median_SALLJ_max_wind_onlySALLJ
y_var = mean_bulk_shear_2_6km_onlySALLJ

filename_x = 'median_SALLJ_max_wind_onlySALLJ'
filename_y = 'mean_bulk_shear_2_6km_onlySALLJ'

xlabel = 'median_SALLJ_max_wind (kts)'
ylabel = 'Mean 2-6 km shear (kts)'

ax.set_aspect(0.5)

############################## Color by SALLJ coverage ###############################

plotted_fig = ax.scatter(x_var,y_var, c=MCS_prop_area_SALLJ_onlySALLJ, cmap='Reds', zorder=2) # MCS_prop_area_SALLJ_onlySALLJ

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

#ax.scatter(x_var, y_var, marker='o', s=90, c=median_SALLJ_max_wind, cmap='Blues', zorder=1)

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

plt.savefig('2%s_vs_%s_scatter%s_%s.png' %(filename_y, filename_x, filter_label, MCS_init_area), dpi=200)

print('saved')
    