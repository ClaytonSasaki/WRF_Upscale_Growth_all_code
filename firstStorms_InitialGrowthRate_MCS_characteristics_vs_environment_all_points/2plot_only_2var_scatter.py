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
import matplotlib.patches as mpatches
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

add_SALLJ = False
add_statistics = False

events_removed = True

plot_type = '2var_scatter_grouped'

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
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/3MCS_characteristics_vs_environment_all_points/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# var choices: q_850, q_flux_850_v, srh_0_1km, srh_0_3km, srh_0_1km_noLatCorrection, srh_0_3km_noLatCorrection, MUCAPE, MUCIN, LCL, LFC, bulk_shear_0_3km, bulk_shear_0_6km, bulk_shear_2_6km

var_plot_str_x = 'q_850' 
var_plot_str_y = 'MUCAPE'

# label choices: 'mean q at 850hPa', 'mean q v flux at 850hPa', '0-1km SRH (m2/s2)', 'MUCAPE (J/Kg)', 'MUCIN (J/kg)', 'LCL (m)', 'LFC (m)', '0-3km Bulk Shear (kts)'
xlabel = 'mean q at 850hPa'
ylabel = 'MUCAPE (J/Kg)'

# Note: to do LFC-LCL, var_plot_str = 'LFC' and uncomment lines labeled 'For LFC-LCL'

stat_type = 'mean' # mean, meadian

# read in files
MCS_prop_area_SALLJ_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
##MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
#MCS_duration_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_duration2_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_majoraxislength_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_growth_stage_time_length_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_growth_stage_time_length_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_ccs_area_growth_filtered_all_points_byEvent_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
median_SALLJ_max_wind_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#median_SALLJ_height_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
var_plot_x_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_%s_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, var_plot_str_x, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
var_plot_y_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_%s_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, var_plot_str_y, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#LCL_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_LCL_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb")) # For LFC-LCL
#MCS_refl_coverage_grt0_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt0_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_refl_coverage_grt30_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt30_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

######## remove certain problematic events ########

if events_removed == True:

    events_bool = np.full(len(var_plot_x_all_points_byEvent), True)

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
    
    events_bool = np.full(len(var_plot_x_all_points_byEvent), True)
    
    events_removed_label = ''
    

#MCS_prop_area_SALLJ_all_points_byEvent = np.array(MCS_prop_area_SALLJ_all_points_byEvent)[events_bool]
#median_SALLJ_max_wind_all_points_byEvent = np.array(median_SALLJ_max_wind_all_points_byEvent)[events_bool]
MCS_ccs_area_growth_filtered_all_points_byEvent = np.array(MCS_ccs_area_growth_filtered_all_points_byEvent)[events_bool]
MCS_growth_stage_time_length_filtered_all_points_byEvent = np.array(MCS_growth_stage_time_length_filtered_all_points_byEvent)[events_bool]
var_plot_x_all_points_byEvent = np.array(var_plot_x_all_points_byEvent)[events_bool]
var_plot_y_all_points_byEvent = np.array(var_plot_y_all_points_byEvent)[events_bool]

#var_plot_all_points_byEvent = var_plot_all_points_byEvent - LCL_all_points_byEvent # For LFC-LCL

# calculate mean/median in variables and save them into 1D array (for MCS and SALLJ characteristics, first value is taken as does not vary spatially)

#MCS_prop_area_SALLJ = np.array([arr[0,0] for arr in MCS_prop_area_SALLJ_all_points_byEvent])
#median_SALLJ_max_wind = np.array([arr[0,0] for arr in median_SALLJ_max_wind_all_points_byEvent])
MCS_ccs_area_growth_filtered = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent])
MCS_growth_stage_time_length_filtered = np.array([arr[0,0] for arr in MCS_growth_stage_time_length_filtered_all_points_byEvent])

if stat_type == 'mean':
    var_plot_x = np.array([np.nanmean(arr) for arr in var_plot_x_all_points_byEvent])
    var_plot_y = np.array([np.nanmean(arr) for arr in var_plot_y_all_points_byEvent])
elif stat_type == 'median':
    var_plot_x = np.array([np.nanmedian(arr) for arr in var_plot_x_all_points_byEvent])
    var_plot_y = np.array([np.nanmedian(arr) for arr in var_plot_y_all_points_byEvent])
else:
    print('Please choose a valid stat_type')

# calculate growth rate

MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
print(MCS_ccs_area_growth_rate_filtered)

# remove events where growth rate is nan (because MCS was found at time 0) or either chosen var is nan

nan_indices = np.isnan(MCS_ccs_area_growth_rate_filtered)

#MCS_prop_area_SALLJ = np.delete(MCS_prop_area_SALLJ, nan_indices)
#median_SALLJ_max_wind = np.delete(median_SALLJ_max_wind, nan_indices)
MCS_ccs_area_growth_filtered = np.delete(MCS_ccs_area_growth_filtered, nan_indices)
MCS_growth_stage_time_length_filtered = np.delete(MCS_growth_stage_time_length_filtered, nan_indices)
MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, nan_indices)
var_plot_x = np.delete(var_plot_x, nan_indices)
var_plot_y = np.delete(var_plot_y, nan_indices)

nan_indices = np.isnan(var_plot_x)

MCS_ccs_area_growth_filtered = np.delete(MCS_ccs_area_growth_filtered, nan_indices)
MCS_growth_stage_time_length_filtered = np.delete(MCS_growth_stage_time_length_filtered, nan_indices)
MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, nan_indices)
var_plot_x = np.delete(var_plot_x, nan_indices)
var_plot_y = np.delete(var_plot_y, nan_indices)

nan_indices = np.isnan(var_plot_y)

MCS_ccs_area_growth_filtered = np.delete(MCS_ccs_area_growth_filtered, nan_indices)
MCS_growth_stage_time_length_filtered = np.delete(MCS_growth_stage_time_length_filtered, nan_indices)
MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, nan_indices)
var_plot_x = np.delete(var_plot_x, nan_indices)
var_plot_y = np.delete(var_plot_y, nan_indices)


# remove incorrect values
if var_plot_str_x == 'srh_0_1km_noLatCorrection':

    incorrect_indices = np.where(var_plot_x < -300)
    var_plot_x = np.delete(var_plot_x, incorrect_indices)
    var_plot_y = np.delete(var_plot_y, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    
elif var_plot_str_y == 'srh_0_1km_noLatCorrection':

    incorrect_indices = np.where(var_plot_y < -300)
    var_plot_x = np.delete(var_plot_x, incorrect_indices)
    var_plot_y = np.delete(var_plot_y, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    
elif var_plot_str_x == 'srh_0_1km':

    incorrect_indices = np.where(var_plot_x < -400)
    var_plot_x = np.delete(var_plot_x, incorrect_indices)
    var_plot_y = np.delete(var_plot_y, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    
elif var_plot_str_y == 'srh_0_1km':

    incorrect_indices = np.where(var_plot_y < -400)
    var_plot_x = np.delete(var_plot_x, incorrect_indices)
    var_plot_y = np.delete(var_plot_y, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    
elif var_plot_str_x == 'LFC':

    incorrect_indices = np.where(var_plot_x > 15000)
    var_plot_x = np.delete(var_plot_x, incorrect_indices)
    var_plot_y = np.delete(var_plot_y, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    
elif var_plot_str_y == 'LFC':

    incorrect_indices = np.where(var_plot_y > 15000)
    var_plot_x = np.delete(var_plot_x, incorrect_indices)
    var_plot_y = np.delete(var_plot_y, incorrect_indices)

    MCS_ccs_area_growth_rate_filtered = np.delete(MCS_ccs_area_growth_rate_filtered, incorrect_indices)
    
else:
    pass


fig, ax = plt.subplots()


if plot_type == 'scatter':

    #------ MCS Length Growth vs 850hPa q -------

    #x_var = mean_q_850
    #y_var = MCS_majoraxislength_growth_filtered
    #
    #filename_x = 'mean_q_850'
    #filename_y = 'MCS_majoraxislength_growth_filtered'
    #
    #xlabel = 'mean q at 850hPa'
    #ylabel = 'MCS_majoraxislength_growth (km/h)'

    #------ MCS Area Growth vs 850hPa q -------

    #x_var = mean_q_850
    #y_var = MCS_ccs_area_growth_filtered
    #
    #filename_x = 'mean_q_850'
    #filename_y = 'MCS_ccs_area_growth_filtered'
    #
    #xlabel = 'mean q at 850hPa'
    #ylabel = 'MCS_ccs_area_growth_filtered (km2/h)'

    #------ MCS Area Growth Rate vs 850hPa q -------

#    x_var = mean_q_850
#    y_var = MCS_ccs_area_growth_rate_filtered
#
#    filename_x = 'mean_q_850'
#    filename_y = 'MCS_ccs_area_growth_rate_filtered'
#
#    xlabel = 'mean q at 850hPa'
#    ylabel = 'ccs growth rate (km2/hr)'

    #------ MCS Area Growth Rate vs 850hPa q flux v -------

    x_var = mean_q_flux_850_v
    y_var = MCS_ccs_area_growth_rate_filtered

    filename_x = 'mean_q_flux_850_v'
    filename_y = 'MCS_ccs_area_growth_rate_filtered'

    xlabel = 'mean q v flux at 850hPa'
    ylabel = 'ccs growth rate (km2/hr)'

    #------ MCS Area Growth Rate vs MUCAPE -------

#    x_var = median_MUCAPE
#    y_var = MCS_ccs_area_growth_rate_filtered
#    
#    filename_x = 'median_MUCAPE'
#    filename_y = 'MCS_ccs_area_growth_rate_filtered'
#    
#    xlabel = 'median MUCAPE'
#    ylabel = 'ccs growth rate (km2/hr)'

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

    #------ MCS Area Growth Rate vs Shear Strength -------

#    x_var = median_bulk_shear_2_6km
#    y_var = MCS_ccs_area_growth_rate_filtered
#    
#    filename_x = 'median_bulk_shear_2_6km'
#    filename_y = 'MCS_ccs_area_growth_rate_filtered'
#    
#    xlabel = 'median 2-6 km shear (kts)'
#    ylabel = 'ccs growth rate (km2/hr)'

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

    # ----------------- ADD SALLJ --------------------

    if add_SALLJ == True:

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

        SALLJ_added_label = '_SALLJ_added'

    else:

        ax.scatter(x_var, y_var, marker='o', s=90)

        SALLJ_added_label = ''

    # ----------------- ADD STATISTICS --------------------

    if add_statistics == True:

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

        stats_added_label = '_stats_added'

    else:

        stats_added_label = ''

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    
    filetype_label = '%s_%s_vs_%s_scatter' %(stat_type, filename_y, filename_x)
    
    flags = '%s%s%s%s' %(filter_label, events_removed_label, SALLJ_added_label, stats_added_label)
    
elif plot_type == '2var_scatter':
    
    plotted_fig = ax.scatter(var_plot_x,var_plot_y, c=MCS_ccs_area_growth_rate_filtered, cmap='Reds', zorder=2)
        
#    slow_growth = MCS_ccs_area_growth_rate_filtered < 5000
#    var_plot_x_slow_growth = var_plot_x[slow_growth]
#    var_plot_y_slow_growth = var_plot_y[slow_growth]
#    
#    ax.scatter(var_plot_x_slow_growth, var_plot_y_slow_growth, marker='o', s=90, color='blue', zorder=1)

    cbar = fig.colorbar(plotted_fig)
    cbar.ax.set_ylabel('ccs growth rate (km2/hr)')
    
    filetype_label = '%s_%s_vs_%s_scatter' %(stat_type, var_plot_str_y, var_plot_str_x)
    
    SALLJ_added_label = ''
    stats_added_label = ''
    
    flags = '%s%s%s%s' %(filter_label, events_removed_label, SALLJ_added_label, stats_added_label)
    
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    
elif plot_type == '2var_scatter_grouped':
    
    plotted_fig = ax.scatter(var_plot_x,var_plot_y, c=MCS_ccs_area_growth_rate_filtered, cmap='Reds', zorder=2)
    
    area_growth_quartile1, area_growth_median, area_growth_quartile3 = np.percentile(MCS_ccs_area_growth_rate_filtered, [25, 50, 75])
    
    slow_growth = MCS_ccs_area_growth_rate_filtered < area_growth_quartile1
    var_plot_x_slow_growth = var_plot_x[slow_growth]
    var_plot_y_slow_growth = var_plot_y[slow_growth]
    
    ax.scatter(var_plot_x_slow_growth, var_plot_y_slow_growth, marker='o', s=90, color='green', zorder=1)
    
    lmedium_growth = np.logical_and(MCS_ccs_area_growth_rate_filtered >= area_growth_quartile1, MCS_ccs_area_growth_rate_filtered < area_growth_median)
    var_plot_x_lmedium_growth = var_plot_x[lmedium_growth]
    var_plot_y_lmedium_growth = var_plot_y[lmedium_growth]
    
    ax.scatter(var_plot_x_lmedium_growth, var_plot_y_lmedium_growth, marker='o', s=90, color='blue', zorder=1)
    
    hmedium_growth = np.logical_and(MCS_ccs_area_growth_rate_filtered >= area_growth_median, MCS_ccs_area_growth_rate_filtered < area_growth_quartile3)
    var_plot_x_hmedium_growth = var_plot_x[hmedium_growth]
    var_plot_y_hmedium_growth = var_plot_y[hmedium_growth]
    
    ax.scatter(var_plot_x_hmedium_growth, var_plot_y_hmedium_growth, marker='o', s=90, color='yellow', zorder=1)
    
    rapid_growth = MCS_ccs_area_growth_rate_filtered >= area_growth_quartile3
    var_plot_x_rapid_growth = var_plot_x[rapid_growth]
    var_plot_y_rapid_growth = var_plot_y[rapid_growth]
    
    ax.scatter(var_plot_x_rapid_growth, var_plot_y_rapid_growth, marker='o', s=90, color='red', zorder=1)

    cbar = fig.colorbar(plotted_fig)
    cbar.ax.set_ylabel('ccs growth rate (km2/hr)')
    
    filetype_label = '%s_%s_vs_%s_scatter_grouped' %(stat_type, var_plot_str_y, var_plot_str_x)
    
    SALLJ_added_label = ''
    stats_added_label = ''
    
    flags = '%s%s%s%s' %(filter_label, events_removed_label, SALLJ_added_label, stats_added_label)
    
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    
    


plt.tight_layout()

print('saving')

specific_outpath = '%sarea_%s%s%s/plots/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

plt.savefig(general_path + specific_outpath + '15%s%s.png' %(filetype_label, flags), dpi=200)

print('saved')
    