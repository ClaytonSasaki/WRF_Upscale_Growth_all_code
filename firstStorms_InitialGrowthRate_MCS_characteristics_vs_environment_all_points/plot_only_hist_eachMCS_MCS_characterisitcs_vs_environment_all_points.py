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
import matplotlib.patches as mpatches
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

plot_rapid_growth_dist = False
plot_slow_growth_dist = False
plot_all_growth_dist = True

remove_refl_points = False
refl_thres = 30

plot_type = 'hist' # violin or hist

####################################################################

# get corresponding file labels usinf chosen inputs
SALLJ_search_text = '__SALLJarea_%s' %(SALLJ_search_area)
env_search_text = '__EnvArea_%s' %(env_search_area)

if MCS_type == 'all_MCS':
    MCS_file_label = 'MCS'
elif MCS_type == 'robustMCS':
    MCS_file_label = 'robustMCS'
else:
    print('MUST choose valid MCS_type')
    
######## Read in data created from processing script ###############

# get input files paths
general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/3MCS_characteristics_vs_environment_all_points/'

specific_inpath = '%sarea_%s%s%s/data/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

# read in files
#MCS_prop_area_SALLJ_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_prop_area_SALLJ_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
##MCS_duration_filtered = pickle.load(open("MCS_duration%s.dat" %(filter_label), "rb")) # NOTE: this duration from 'mcs_length' variable is all zeros
#MCS_duration_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_duration2_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_majoraxislength_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_majoraxislength_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_ccs_area_growth_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_growth_stage_time_length_filtered_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_growth_stage_time_length_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_ccs_area_growth_filtered_all_points_byEvent_2hr = pickle.load(open(general_path + specific_inpath + "%s_ccs_area_growth_2hr_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_3km_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_3km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_0_6km_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_bulk_shear_0_6km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
bulk_shear_2_6km_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_bulk_shear_2_6km_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#median_SALLJ_max_wind_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_median_SALLJ_max_wind_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#median_SALLJ_height_all_points_byEvent = pickle.load(open(general_path + specific_inpath + "%s_median_SALLJ_height_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_q_850_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_q_850_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
MCS_MUCAPE_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_MUCAPE_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_wv_flux_850_v_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_wv_flux_850_v_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_refl_coverage_grt0_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt0_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))
#MCS_refl_coverage_grt30_byEvent = pickle.load( open(general_path + specific_inpath + "%s_refl_coverage_grt30_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

for event_num in range(len(bulk_shear_0_3km_all_points_byEvent)):

    # get array for each event seperately

    #MCS_prop_area_SALLJ = MCS_prop_area_SALLJ_all_points_byEvent[i]
    #MCS_duration_filtered = MCS_duration_filtered_all_points_byEvent[i]
    #MCS_majoraxislength_growth_filtered = MCS_majoraxislength_growth_filtered_all_points_byEvent[i]
    #MCS_ccs_area_growth_filtered_spreadByEvent = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent])
    #MCS_ccs_area_growth_filtered_spreadByEvent_2hr = np.array([arr[0,0] for arr in MCS_ccs_area_growth_filtered_all_points_byEvent_2hr])
    #print('MCS_ccs_area_growth_filtered_spreadByEvent', MCS_ccs_area_growth_filtered_spreadByEvent)
    bulk_shear_0_3km = bulk_shear_0_3km_all_points_byEvent[event_num]
    bulk_shear_0_6km = bulk_shear_0_6km_all_points_byEvent[event_num]
    bulk_shear_2_6km = bulk_shear_2_6km_all_points_byEvent[event_num]
    MCS_q_850 = MCS_q_850_all_points_byEvent[event_num]
    MCS_MUCAPE = MCS_MUCAPE_all_points_byEvent[event_num]
    #median_SALLJ_max_wind = median_SALLJ_max_wind_all_points_byEvent[i]
    #median_SALLJ_height = median_SALLJ_height_all_points_byEvent[i]

    #q_850_quart75ByEvent = np.array([np.percentile(arr, 75) for arr in MCS_q_850_all_points_byEvent])
    #q_850_quart90ByEvent = np.array([np.percentile(arr, 90) for arr in MCS_q_850_all_points_byEvent])
    #MCS_MUCAPE_quart90ByEvent = np.array([np.percentile(arr, 90) for arr in MCS_MUCAPE_all_points_byEvent])
    #MCS_wv_flux_quart90ByEvent = np.array([np.percentile(arr, 90) for arr in MCS_wv_flux_850_v_all_points_byEvent])
    #
    #q_850_mean = np.array([np.nanmean(arr) for arr in MCS_q_850_all_points_byEvent])
    #MCS_MUCAPE_mean = np.array([np.nanmean(arr) for arr in MCS_MUCAPE_all_points_byEvent])
    #MCS_wv_flux_850_mean = np.array([np.nanmean(arr) for arr in MCS_wv_flux_850_v_all_points_byEvent])
    MCS_ccs_area_growth_rate = MCS_ccs_area_growth_filtered_all_points_byEvent[event_num][0,0]/MCS_growth_stage_time_length_filtered_all_points_byEvent[event_num][0,0]

    #### used to remove points w/reflectivity over certain threshold ####

    MCS_max_refl_all_points_byEvent = pickle.load( open(general_path + specific_inpath + "%s_max_refl_all_points_byEvent%s_%s_%s_%s.dat" %(MCS_file_label, filter_label, MCS_init_area, SALLJ_search_area, env_search_area), "rb"))

    #print(bulk_shear_0_3km.size)
    #print(MCS_max_refl_all_points_byEvent.size)

    ###### remove points with reflectivity greater than threshold #######

    if remove_refl_points == True:

        # filter each environmental var to only be points with refl < refl_thres
        bulk_shear_0_3km = bulk_shear_0_3km[MCS_max_refl_all_points_byEvent < refl_thres]
        bulk_shear_0_6km = bulk_shear_0_6km[MCS_max_refl_all_points_byEvent < refl_thres]
        bulk_shear_2_6km = bulk_shear_2_6km[MCS_max_refl_all_points_byEvent < refl_thres]
        MCS_q_850 = MCS_q_850[MCS_max_refl_all_points_byEvent < refl_thres]
        MCS_MUCAPE = MCS_MUCAPE[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_wv_flux_850_v = MCS_wv_flux_850_v[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_q_flux_850 = MCS_q_flux_850[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_q_flux_850_v = MCS_q_flux_850_v[MCS_max_refl_all_points_byEvent < refl_thres]

    #    # filter MCS characteristic var to match size of env var
    #    MCS_prop_area_SALLJ = MCS_prop_area_SALLJ[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_duration_filtered = MCS_duration_filtered[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_majoraxislength_growth_filtered = MCS_majoraxislength_growth_filtered[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_ccs_area_growth_filtered = MCS_ccs_area_growth_filtered[MCS_max_refl_all_points_byEvent < refl_thres]
    #    MCS_growth_stage_time_length_filtered = MCS_growth_stage_time_length_filtered[MCS_max_refl_all_points_byEvent < refl_thres]
    #    #MCS_ccs_area_growth_filtered_2hr = MCS_ccs_area_growth_filtered_2hr[MCS_max_refl_all_points < refl_thres]
    #    median_SALLJ_max_wind = median_SALLJ_max_wind[MCS_max_refl_all_points_byEvent < refl_thres]
    #    median_SALLJ_height = median_SALLJ_height[MCS_max_refl_all_points_byEvent < refl_thres]

        MCS_max_refl_all_points_byEvent = MCS_max_refl_all_points_byEvent[MCS_max_refl_all_points_byEvent < refl_thres] # filter max_refl to match others but only after done using it as filter

        dBZ_label = '_removed_grt_30_dBZ'

    else:

        dBZ_label = ''

    #MCS_ccs_area_growth_rate_filtered = MCS_ccs_area_growth_filtered/MCS_growth_stage_time_length_filtered
    #median_MCS_ccs_area_growth_rate_filtered = np.nanmedian(MCS_ccs_area_growth_rate_filtered)
    #q75_MCS_ccs_area_growth_rate_filtered = np.nanpercentile(MCS_ccs_area_growth_rate_filtered, 75)
    #q25_MCS_ccs_area_growth_rate_filtered = np.nanpercentile(MCS_ccs_area_growth_rate_filtered, 25)

    ############ Change variable plotting ############ 

    data = bulk_shear_2_6km
    variable_name = 'bulk_shear_2_6km' # MCS_prop_area_SALLJ, bulk_shear_0_3km, MCS_q_850, MCS_MUCAPE, MCS_wv_flux_850_v ...
    var_label = 'bulk 2-6 km shear' # proportion of area w/SALLJ, bulk 0-3 km shear, 850-hPa q, MUCAPE (J/kg), 850hPa v wv flux (m/s) ...

    ############# Change variable for seperation ############
    #
    #var_seperation = MCS_ccs_area_growth_rate_filtered
    #var_seperation_name = 'ccs_area_growth_rate'
    #seperation_threshold = median_MCS_ccs_area_growth_rate_filtered

    ############# Seperate variable for plotted by seperation variable ##############

    # remove nans because violin plot cannot handle nans
    data_nonan = data[~np.isnan(data)]
    #var_seperation = var_seperation[~np.isnan(data)]

    #mask_rapid_growth = var_seperation >= seperation_threshold
    #mask_slow_growth = var_seperation < seperation_threshold

    data_all_growth = data_nonan
    #data_rapid_growth = data_nonan[mask_rapid_growth]
    #data_slow_growth = data_nonan[mask_slow_growth]

    ########################## OLD VERISON OF GETTING VAR AND FILTERING #################################
    #mask1 = var_seperation >= seperation_threshold
    #
    #MCS_prop_area_SALLJ_mask1 = MCS_prop_area_SALLJ[mask1]
    #MCS_duration_filtered_mask1 = MCS_duration_filtered[mask1]
    #MCS_majoraxislength_growth_filtered_mask1 = MCS_majoraxislength_growth_filtered[mask1]
    #MCS_ccs_area_growth_filtered_mask1 = MCS_ccs_area_growth_filtered[mask1]
    #MCS_ccs_area_growth_filtered_2hr_mask1 = MCS_ccs_area_growth_filtered_2hr[mask1]
    #bulk_shear_0_3km_mask1 = bulk_shear_0_3km[mask1]
    #bulk_shear_0_6km_mask1 = bulk_shear_0_6km[mask1]
    #bulk_shear_2_6km_mask1 = bulk_shear_2_6km[mask1]
    #median_SALLJ_max_wind_mask1 = median_SALLJ_max_wind[mask1]
    #median_SALLJ_height_mask1 = median_SALLJ_height[mask1]
    #MCS_q_850_mask1 = MCS_q_850[mask1]
    #MCS_MUCAPE_mask1 = MCS_MUCAPE[mask1]
    ##MCS_wv_flux_850_mask1 = MCS_wv_flux_850[mask1]
    #MCS_wv_flux_850_v_mask1 = MCS_wv_flux_850_v[mask1]
    #MCS_q_flux_850_mask1 = MCS_q_flux_850[mask1]
    #MCS_q_flux_850_v_mask1 = MCS_q_flux_850_v[mask1]
    #
    #print('# rapid growth MCS', len(MCS_ccs_area_growth_filtered_2hr_mask1))
    #
    #mask2 = var_seperation < seperation_threshold
    #
    #MCS_prop_area_SALLJ_mask2 = MCS_prop_area_SALLJ[mask2]
    #MCS_duration_filtered_mask2 = MCS_duration_filtered[mask2]
    #MCS_majoraxislength_growth_filtered_mask2 = MCS_majoraxislength_growth_filtered[mask2]
    #MCS_ccs_area_growth_filtered_mask2 = MCS_ccs_area_growth_filtered[mask2]
    #MCS_ccs_area_growth_filtered_2hr_mask2 = MCS_ccs_area_growth_filtered_2hr[mask2]
    #bulk_shear_0_3km_mask2 = bulk_shear_0_3km[mask2]
    #bulk_shear_0_6km_mask2 = bulk_shear_0_6km[mask2]
    #bulk_shear_2_6km_mask2 = bulk_shear_2_6km[mask2]
    #median_SALLJ_max_wind_mask2 = median_SALLJ_max_wind[mask2]
    #median_SALLJ_height_mask2 = median_SALLJ_height[mask2]
    #MCS_q_850_mask2 = MCS_q_850[mask2]
    #MCS_MUCAPE_mask2 = MCS_MUCAPE[mask2]
    ##MCS_wv_flux_850_mask2 = MCS_wv_flux_850[mask2]
    #MCS_wv_flux_850_v_mask2 = MCS_wv_flux_850_v[mask2]
    #MCS_q_flux_850_mask2 = MCS_q_flux_850[mask2]
    #MCS_q_flux_850_v_mask2 = MCS_q_flux_850_v[mask2]
    #
    #print('# slow growth MCS', len(MCS_ccs_area_growth_filtered_2hr_mask2))

    ############ Change variable plotting ############ 

    #data_all_growth = bulk_shear_0_3km
    #data_rapid_growth = bulk_shear_0_3km_mask1
    #data_slow_growth = bulk_shear_0_3km_mask2
    #variable_name = 'bulk_shear_0_3km' # MCS_prop_area_SALLJ, bulk_shear_0_3km, MCS_q_850, MCS_MUCAPE, MCS_wv_flux_850_v ...
    #var_label = 'bulk 0-3 km shear' # proportion of area w/SALLJ, bulk 0-3 km shear, 850-hPa q, MUCAPE (J/kg), 850hPa v wv flux (m/s) ...

    #################################################################################################### 

    ########################## Plotting ##############################
    fig, ax = plt.subplots()
    labels = []
    labels_text_only = []

    data_list = []
    color_list = []
    
    if ~np.isnan(MCS_ccs_area_growth_rate):
    
        ax.text(0.01, 0.01, '%d km/hr' %(MCS_ccs_area_growth_rate), transform=ax.transAxes)
        
    else:
        
        ax.text(0.01, 0.01, 'NaN', transform=ax.transAxes)

    # setting data and labels by plots chosen
    if plot_all_growth_dist == True:

        data_list.append(data_all_growth)
        labels.append((mpatches.Patch(color='gray', alpha=0.2), 'all MCS'))
        labels_text_only.append('all MCS')

        color_list.append('gray')
        dist_text_all_growth = '_all' # set text for filename

    else:

        dist_text_all_growth = ''

    if plot_rapid_growth_dist == True:

        data_list.append(data_rapid_growth)
        labels.append((mpatches.Patch(color='blue', alpha=0.2), 'rapid growth MCS'))
        labels_text_only.append('rapid growth MCS')

        color_list.append('blue')
        dist_text_rapid_growth = '_rapid' # set text for filename

    else:

        dist_text_rapid_growth = ''

    if plot_slow_growth_dist == True:

        data_list.append(data_slow_growth)
        labels.append((mpatches.Patch(color='red', alpha=0.2), 'slow growth MCS'))
        labels_text_only.append('slow growth MCS')

        color_list.append('red')
        dist_text_slow_growth = '_slow' # set text for filename

    else:

        dist_text_slow_growth = ''

    # plotting based on plot type chosen  
    if plot_type == 'hist':

        for i, data in enumerate(data_list):

            median_value = np.nanmedian(data)  
            plt.hist(data, bins=30, alpha=0.3, density=True, label=labels_text_only[i], color=color_list[i])
            plt.axvline(median_value, color=color_list[i], linestyle='dashed', linewidth=2, label=f'Median: {median_value:.2f}')

            # adding labels
            plt.xlabel(var_label)
            plt.ylabel('Probability Density')

            # show legend
            plt.legend()

    elif plot_type == 'violin':

        alpha = 0.2

        def adjacent_values(vals, q1, q3):
            upper_adjacent_value = q3 + (q3 - q1) * 1.5
            upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

            lower_adjacent_value = q1 - (q3 - q1) * 1.5
            lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
            return lower_adjacent_value, upper_adjacent_value

        for i, data in enumerate(data_list):

            violin = plt.violinplot(data, showmeans=False, showmedians=False, showextrema=False)

            for pc in violin['bodies']:
                pc.set_facecolor(color_list[i])
                #pc.set_edgecolor('black')
                pc.set_alpha(alpha)

            #violin['cmedians'].set_color(color_list[i])
            #violin['cmedians'].set_linestyle('--') # Set the line style to dashed

            # shows quartiles and median, comment out lines above (violin['cmedians']) and set 'showmedians' above to False to use 
            quartile1, median, quartile3 = np.percentile(data, [25, 50, 75])
            whiskers = np.array([adjacent_values(data, quartile1, quartile3)])
            whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]

            ax.scatter(1, median, marker='o', color=color_list[i], s=50, zorder=3)
            ax.vlines(1, quartile1, quartile3, color=color_list[i], linestyle='-', lw=5, alpha=0.8)
            #ax.vlines(1, whiskers_min, whiskers_max, color='k', linestyle='-', lw=1)

            # adding labels
            plt.ylabel(var_label)

            # adding legend
            plt.legend(*zip(*labels), loc=2)

    plt.tight_layout()

    print('saving')

    specific_outpath = '%sarea_%s%s%s/plots_each_event/' %(MCS_file_label, MCS_init_area, SALLJ_search_text, env_search_text)

    plt.savefig(general_path + specific_outpath + '%d_%s_%s_%s%s.png' %(event_num, variable_name, plot_type, filter_label, dBZ_label), dpi=200)

    print('saved')
    
    plt.close()
    