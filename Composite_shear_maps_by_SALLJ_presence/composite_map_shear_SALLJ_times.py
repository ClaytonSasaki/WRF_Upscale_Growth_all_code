#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Sep 13 2024

@author: crs326

Creates .dat files and plots of composite shear for all SALLJ times

Based upon 'relaxed_SALLJ_spatial_map_frequency_efficient_v4.py', 'plot_only_SALLJ_spatial_map_frequency_efficient_v4.py', 'plot_only_relaxed_SALLJ_distributions_by_height.py' among other files

Uses '3relaxed_SALLJ_times_VMRS_1101_0430.dat' file copied over. See README for more.

"""

import matplotlib
matplotlib.use('Agg') 

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import datetime
import pickle
import math
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import pandas as pd
import os
from scipy import stats
from scipy.interpolate import interp1d
import xarray
#import time

# ------------------------- cmaps ------------------------------

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

import matplotlib.colors as colors

#start_time = time.time()

# used for contourf of terrain 
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap = plt.get_cmap('terrain')
new_cmap = truncate_colormap(cmap, 0.21, 1)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

# used for pcolormesh of wind
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('gray')
new_cmap2 = truncate_colormap(cmap2, 0.3, 0.75)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# used for contour of wind
#def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
#    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
#    return new_cmap
#
#arr = np.linspace(0, 50, 100).reshape((10, 10))
#fig, ax = plt.subplots(ncols=2)
#
#cmap2 = plt.get_cmap('gnuplot2')
#new_cmap2 = truncate_colormap(cmap2, 0.2, 1.0)
#ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
#ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# ----------------- reading in file/plotting terrain -----------------

outpath = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/'

# Open the NetCDF file
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

#change the start and end days to be plotted
start_date = 20181101
end_date = 20190430

########### for finding LLJs ###########

# 2: high - allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow
crit = [2 ,3200, 5700, 19.4384, 7.77538] # relaxed to 10 m/s max and 4 m/s decrease

# min and max pressure to plot
min_pres = 450
max_pres = 975

crit_num = crit[0]
max_search_hgt = crit[1]
min_search_hgt = crit[2]
max_wind_threshold = crit[3]
decrease_to_min_threshold = crit[4]

######################################

# get dates from times to only run for correct day folders
input_start_day = pd.to_datetime(str(start_date), format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(end_date), format='%Y%m%d', errors='ignore')

SALLJ_location = 'VMRS'

Calculate = False
Plot = True

if Calculate == True:

    SALLJ_time_list = pickle.load(open("3relaxed_SALLJ_times_%s_%s_%s.dat" %(SALLJ_location, start_date, end_date), "rb"))

    SALLJ_time_list_arr = np.array([str(SALLJ_time) for SALLJ_time in SALLJ_time_list])

    SALLJ_time_list_arr_dt = [datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in SALLJ_time_list_arr]

    variable_spatial_SALLJ_times = np.empty((1356, 499, 599)) # 1356 should be replaced with len(SALLJ_time_list)
    variable_spatial_SALLJ_times[:] = np.NaN

    SALLJ_count = 0

    plot_terrain = False

    for SALLJ_time in SALLJ_time_list_arr_dt:

        dt_date_str = SALLJ_time.date().strftime('%Y%m%d')

        file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(SALLJ_time.strftime('%Y'), SALLJ_time.strftime('%m'), SALLJ_time.strftime('%d'), SALLJ_time.strftime('%H'))

        file_path = path_wrf + dt_date_str + file_name

        print(file_path)

        # get the netCDF
        ncfile = Dataset(file_path,'r')

        # plot terrain only once
        if(plot_terrain == True):

            # Get the terrain height
            terrain = getvar(ncfile, 'HGT')

            # Get the latitude and longitude points
            lats, lons = latlon_coords(terrain)

            # Get the cartopy mapping object
            cart_proj = get_cartopy(terrain)

            # save terrain, lats, lons, cart_proj 
            pickle.dump(terrain, open("terrain.dat", "wb"))

            pickle.dump(lats, open("lats.dat", "wb"))

            pickle.dump(lons, open("lons.dat", "wb"))

            pickle.dump(cart_proj, open("cart_proj.dat", "wb"))

            # Create a figure
            fig = plt.figure(figsize=(12,6))
            # Set the GeoAxes to the projection used by WRF
            ax = plt.axes(projection=cart_proj)

            #ax.add_feature(LAND)
            ax.add_feature(OCEAN, zorder=2)
            #ax.add_feature(COASTLINE)
            ax.add_feature(LAKES, alpha=0.5, zorder=2)
            ax.add_feature(RIVERS, zorder=2)
            ax.add_feature(BORDERS, edgecolor='gray', zorder=2)

            # Make  filled contours for the terrain
            terrain_plot = plt.contourf(to_np(lons), to_np(lats), to_np(terrain), levels=[0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000,5250,5500,5750,6000,6250], extend='max', transform=crs.PlateCarree(), cmap=new_cmap, zorder=1)

            # Add a color bar
            cbar = plt.colorbar(terrain_plot, ax=ax, shrink=.98)
            cbar.set_label('Elevation (m)')

            # Set the map bound

        #    area = 'close_zoom'
        #    lat_bottom_left = -30.0
        #    lon_bottom_left = -64.0
        #    
        #    lat_top_right = -29.0
        #    lon_top_right = -63.0

        #    area = 'medium_zoom'
        #    lat_bottom_left = -32.5
        #    lon_bottom_left = -66.0
        #    
        #    lat_top_right = -27.5
        #    lon_top_right = -60.0

        #    area = 'less_zoom'
        #    lat_bottom_left = -36.0
        #    lon_bottom_left = -68.0
        #    
        #    lat_top_right = -21.0
        #    lon_top_right = -56.0

            # when using all area use three lines below and comment out bounds, ax.set_..., and bottom_left_xy and top_right_xy
            area = 'all_area'   
            ax.set_xlim(cartopy_xlim(terrain))
            ax.set_ylim(cartopy_ylim(terrain))

        #    bounds = GeoBounds(CoordPair(lat=lat_bottom_left, lon=lon_bottom_left),
        #                           CoordPair(lat=lat_top_right, lon=lon_top_right))
        #    
        #    ax.set_xlim(cartopy_xlim(terrain, geobounds=bounds))
        #    ax.set_ylim(cartopy_ylim(terrain, geobounds=bounds))
        #
        #    bottom_left_xy = ll_to_xy(ncfile, lat_bottom_left, lon_bottom_left)
        #    top_right_xy = ll_to_xy(ncfile, lat_top_right, lon_top_right) 

            # Add the gridlines
            gl = ax.gridlines(crs=crs.PlateCarree(), linewidth=1, color='gray', alpha=0.5, draw_labels=True, linestyle='--', zorder=5)
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = True
            gl.bottom_labels = True

            ax.scatter(-64.212, -31.298, s=30, color='gray', transform=crs.PlateCarree(), zorder=6)
            ax.scatter(-63.726, -29.906, s=30, color='gray', transform=crs.PlateCarree(), zorder=6)

            # set plot_terrain to false as only need to complete once
            plot_terrain = False

            print('plotted terrain')


        ############################# read in some variables #############################
        print('reading in variables for: %s %s'  %(dt_date_str, SALLJ_time.strftime('%H')))
        # read in the variables 
        pres = getvar(ncfile, 'pressure')
        hght = getvar(ncfile, 'height', msl=False, units='m') # in m NOTE: this is AGL not MSL!!
        u = getvar(ncfile, 'ua', units='kt') # in kts
        v = getvar(ncfile, 'va', units='kt') # in kts
        speed, drct = getvar(ncfile, 'wspd_wdir', units='kt') # in kts

        #################### get winds at chosen heights

        level = 500

        u_level = np.squeeze(interplevel(u, pres, level))
        v_level = np.squeeze(interplevel(v, pres, level))

        u_sfc = np.squeeze(u[3, :, :])
        v_sfc = np.squeeze(v[3, :, :])

        u_diff = u_level - u_sfc
        v_diff = v_level - v_sfc

        mag_diff = np.sqrt(u_diff**2 + v_diff**2)

        #smooth_mag_diff = smooth2d(mag_diff, 10)

        #mask = terrain <= 3000
        #mag_diff_terrainLess3000 = np.where(mask, mag_diff, np.nan)

        variable_spatial_SALLJ_times[SALLJ_count,:,:] = mag_diff

        #mean_variable_spatial_SALLJ_times = np.nanmean(variable_spatial_SALLJ_times, axis=0)

        SALLJ_count = SALLJ_count + 1

        #print(variable_spatial_SALLJ_times[:,0,0])
        #print(mean_variable_spatial_SALLJ_times)

    mean_variable_spatial_SALLJ_times = np.nanmean(variable_spatial_SALLJ_times, axis=0)
    print(mean_variable_spatial_SALLJ_times)

    pickle.dump(mean_variable_spatial_SALLJ_times, open("mean_%s_SALLJ_times_%s_%s.dat" %('0_6km_shear', start_date, end_date), "wb"))

################################# plot #################################

if Plot == True:

    terrain = pickle.load(open("terrain.dat", "rb"))
    lats = pickle.load(open("lats.dat", "rb"))
    lons = pickle.load(open("lons.dat", "rb"))
    cart_proj = pickle.load(open("cart_proj.dat", "rb"))

    # Create a figure
    fig = plt.figure(figsize=(12,6))
    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)

    #ax.add_feature(LAND)
    ax.add_feature(OCEAN, zorder=2)
    #ax.add_feature(COASTLINE)
    #ax.add_feature(LAKES, alpha=0.5, zorder=2)
    #ax.add_feature(RIVERS, zorder=2)
    ax.add_feature(BORDERS, edgecolor='darkgray', zorder=2)

    # Make  filled contours for the terrain
    terrain_plot = plt.contourf(to_np(lons), to_np(lats), to_np(terrain), levels=[0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000], extend='max', transform=crs.PlateCarree(), cmap=new_cmap2, zorder=1)

    # Outline terrain greater than 500 m
    #terrain_grt_500_plot = plt.contour(to_np(lons), to_np(lats), to_np(terrain), levels=[500,750,1000,1250], extend='max', transform=crs.PlateCarree(), cmap=new_cmap2, zorder=5)

    # Add a color bar for terrain
    cbar = plt.colorbar(terrain_plot, ax=ax, shrink=.98)
    cbar.set_label('Elevation (m)')

    area = 'all_area'   
    ax.set_xlim(cartopy_xlim(terrain))
    ax.set_ylim(cartopy_ylim(terrain))

    # Add the gridlines
    gl = ax.gridlines(crs=crs.PlateCarree(), linewidth=1, color='gray', alpha=0.5, draw_labels=True, linestyle='--', zorder=5)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = True
    gl.bottom_labels = True

    ax.scatter(-64.212, -31.298, s=30, color='black', transform=crs.PlateCarree(), zorder=6)
    ax.scatter(-63.726, -29.906, s=30, color='black', transform=crs.PlateCarree(), zorder=6)


    print('plotted terrain')

    mean_variable_spatial_SALLJ_times = pickle.load(open("mean_%s_SALLJ_times_%s_%s.dat" %('0_6km_shear', start_date, end_date), "rb"))

    mean_variable_spatial_SALLJ_times_smooth = smooth2d(mean_variable_spatial_SALLJ_times, 6, cenweight=4)

    levels = np.arange(5,55,5)

    shear_0_6km = plt.contour(lons, lats, mean_variable_spatial_SALLJ_times_smooth, levels=levels, cmap='Reds', transform=crs.PlateCarree())

    ax.clabel(shear_0_6km, shear_0_6km.levels, inline=True, fmt='%.0f')

    plt.tight_layout()

    print('saving')

    plt.savefig(outpath + 'shear_0_6km_composite_map_SALLJ_times_%s_%s_%s.png' %(start_date, end_date, SALLJ_location))

    print('saved')

    
    #shear_0_6km = plt.contour(lons,lats,mag_diff_terrainLess3000, levels=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100], cmap='Reds', transform=crs.PlateCarree())
    
    #ax.clabel(shear_0_6km, shear_0_6km.levels, inline=True, fmt='%.0f')
    
    #plt.savefig(outpath + 'shear_0_6km_map_%s_%s.png' %(dt_date_str, SALLJ_time.strftime('%H')))

#    ############################# find SALLJ #############################
#
#    # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind minimum 
#    interp_levels_min = np.arange(0,min_search_hgt+250,250)
#
#    pres_below_min = interplevel(pres, hght, interp_levels_min)
#    hght_below_min = interplevel(hght, hght, interp_levels_min)
#    speed_below_min = interplevel(speed, hght, interp_levels_min)
#    drct_below_min = interplevel(drct, hght, interp_levels_min)
#
#    # interpolate variables to 250 meter vertical spacing up to the search heigth for the wind maximum
#    interp_levels_max = np.arange(0,max_search_hgt+250,250)
#
#    pres_below_max = interplevel(pres, hght, interp_levels_max)
#    hght_below_max = interplevel(hght, hght, interp_levels_max)
#    speed_below_max = interplevel(speed, hght, interp_levels_max)
#    drct_below_max = interplevel(drct, hght, interp_levels_max)
#
#    print(pres_below_max.shape)
#    print(hght_below_max.shape)
#
#    point_to_chceck_xy = ll_to_xy(ncfile, -29.906, -63.726)
#
#    print('pres_below_max', pres_below_max[:,point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#    #print('hght_below_max', hght_below_max[:,point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#    print('speed_below_max', speed_below_max[:,point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#
#    ################ with xarray ###################
#
#    # get max wind
#    max_wind = speed_below_max.max(dim='level')
#
#    print('max_wind', max_wind[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#
#    # get height of max wind
#    level_max_wind = hght_below_max.isel(level=speed_below_max.argmax('level'))
##                level_max_wind = speed_below_max.idxmax('level')
##                level_max_wind = speed_below_max[np.where(np.equal(speed_below_max,max_wind))].get_index('level')
#
#    #print('level_max_wind', level_max_wind.values)
#
#    # get pressure at max wind
#    pres_max_wind = pres_below_max.isel(level=speed_below_max.argmax('level'))
#
#    # get direction at max wind
#    drct_at_max_wind = drct_below_max.isel(level=speed_below_max.argmax('level'))
#
#    #print('drct_at_max_wind', drct_at_max_wind)
#
#    # set dim of level of max wind array to 'level' (height)
#    level_max_wind_subtract = xarray.DataArray(level_max_wind, dims=['south_north', 'west_east'])
#
#    #print('level_max_wind_subtract', level_max_wind_subtract)
#
#    #print('hght_below_max', hght_below_max.values)
#
#    # subtract the heights of max wind from the heights
#    hght_minus_level_max_wind = hght_below_min - level_max_wind_subtract
#
#    #print('hght_minus_level_max_wind', hght_minus_level_max_wind)
#
#    speed_below_min_masked_below_max = speed_below_min.where(hght_minus_level_max_wind > 0., np.nan) 
#
#    #print('speed_below_min_masked_below_max', speed_below_min_masked_below_max)
#
#    # get min wind above max
#    min_wind = speed_below_min_masked_below_max.min(dim='level')
#
#    print('min_wind', min_wind[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#
#    # checks if max_wind meets threshold and keeps value if it does meet the threshold
#    max_wind_meeting_threshold = max_wind.where(max_wind > max_wind_threshold, np.nan)
#
#    print('max_wind_meeting_threshold', max_wind_meeting_threshold[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#
#    #print('max_wind_meeting_threshold', max_wind_meeting_threshold)
#
#    # calculates decrease to min wind
#    decrease_to_min = max_wind_meeting_threshold - min_wind
#
#    #print('decrease_to_min', decrease_to_min)
#
#    # checks if decrease_to_min meets threshold and keeps value if it does meet the threshold
#    decrease_to_min_meeting_threshold = decrease_to_min.where(decrease_to_min > decrease_to_min_threshold, np.nan)
#
#    #print('decrease_to_min_meeting_threshold', decrease_to_min_meeting_threshold)
#
#    # checks to see if the values met the other criteria (was not replaced with nan) and if it does leave the value
#    drct_at_max_wind_meeting_threshold = drct_at_max_wind.where(np.isnan(decrease_to_min_meeting_threshold) == False, np.nan)
#
#    print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#
#    # checks to see if wind at max_wind is from a northerly direction and keeps value if it is
#    drct_at_max_wind_meeting_threshold = drct_at_max_wind_meeting_threshold.where((drct_at_max_wind_meeting_threshold <= 45) | (drct_at_max_wind_meeting_threshold >= 315), np.nan)
#
#    #print('drct_at_max_wind_meeting_threshold', drct_at_max_wind_meeting_threshold)
#
#    # get pressure of max_wind of points that meet SALLJ criteria
#    pres_max_wind_SALLJs = pres_max_wind.where(np.isnan(drct_at_max_wind_meeting_threshold) == False, np.nan)
#
#    print('pres_max_wind_SALLJs', pres_max_wind_SALLJs[point_to_chceck_xy[1],point_to_chceck_xy[0]].values)
#
#    #pres_max_wind_SALLJs[np.isnan(pres_max_wind_SALLJs) == False] = 1
#
#    pres_max_wind_SALLJs = pres_max_wind_SALLJs.where(np.isnan(drct_at_max_wind_meeting_threshold) == True, 1)
#
#    pres_max_wind_SALLJs_for_all_times = np.nansum(np.dstack((pres_max_wind_SALLJs,pres_max_wind_SALLJs_for_all_times)),2)
#
#    ncfile.close()
#
#    hour_counter = hour_counter + 1
#
#print('pres_max_wind_SALLJs_for_all_times', pres_max_wind_SALLJs_for_all_times[point_to_chceck_xy[1],point_to_chceck_xy[0]])
#
##pickle.dump(pres_max_wind_SALLJs_for_all_times, open("pre_SALLJ_all_times.dat", "wb"))
##
##pres_max_wind_SALLJs_for_all_times[np.isnan(pres_max_wind_SALLJs_for_all_times) == False] = 1
##
##print('pres_max_wind_SALLJs_for_all_times', pres_max_wind_SALLJs_for_all_times[:,point_to_chceck_xy[1],point_to_chceck_xy[0]])
##binary_SALLJ = np.nansum(pres_max_wind_SALLJs_for_all_times, axis=0)
##
###pres_max_wind_SALLJs_for_all_times = pres_max_wind_SALLJs_for_all_times.where(np.isnan(pres_max_wind_SALLJs_for_all_times) == True, 1)
###binary_SALLJ = pres_max_wind_SALLJs_for_all_times.sum(dim=0)
##print('binary_SALLJ', binary_SALLJ[point_to_chceck_xy[1],point_to_chceck_xy[0]])
##print(np.amax(binary_SALLJ))
#
#print(np.amax(pres_max_wind_SALLJs_for_all_times))
#
#pickle.dump(pres_max_wind_SALLJs_for_all_times, open("NEW_relaxed_SALLJ_frequency_v4_%s_%s.dat" %(start_date, end_date), "wb"))
#
## For 2 days (to get levels 0 to 55 by 5%)
##plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[2.4,4.8,7.2,9.6,12,14.4,16.8,19.2,21.6,24,26.4], extend='max', cmap='Greys', zorder=4, transform=crs.PlateCarree())
#
## For 30 days (to get levels 0 to 55 by 5%)
##plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[36,72,108,144,180,216,252,288,324,360,396], extend='max', cmap='Greys', zorder=4, transform=crs.PlateCarree())
#
## For 61 days (to get levels 0 to 55 by 5%)
#plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[217.2,434.4,651.6,868.8,1086,1303.2,1520.4,1737.6,1954.8,2172,2389.2], extend='max', cmap='Greys', zorder=4, transform=crs.PlateCarree())
#
#
### write data to .csv file
##with open('PNNL_WRF_SDC_CAPE_time_series_%s_%s_%s.csv' % (area_names[i], str(start_date)[4:], str(end_date)[4:]), 'w') as f:
##    writer = csv.writer(f)
##
##    # write the times to a row
##    writer.writerow(time_plot)
##
##    # write the data to the next rows
##    writer.writerow(Average_MUCAPE_subset_list)
##
##print('done')
#
#
#
##    # get the pressure of max_wind of points that meet SALLJ criteria that fall into pressure groupings
##    pres_max_wind_SALLJs_less700 = pres_max_wind_SALLJs.where(pres_max_wind_SALLJs <= 700.0, np.nan)
##    pres_max_wind_SALLJs_700_800 = pres_max_wind_SALLJs.where((pres_max_wind_SALLJs > 700.0) & (pres_max_wind_SALLJs <= 800.0), np.nan)
##    pres_max_wind_SALLJs_800_900 = pres_max_wind_SALLJs.where((pres_max_wind_SALLJs > 800.0) & (pres_max_wind_SALLJs <= 900.0), np.nan)
##    pres_max_wind_SALLJs_grt900 = pres_max_wind_SALLJs.where(pres_max_wind_SALLJs >= 900.0, np.nan)
##                            
##                            if max_wind_pres <= 700:
##                                height_index = 'gold'
##                            elif max_wind_pres > 700 and max_wind_pres <=750:
##                                height_index = 'orange'
##                            elif max_wind_pres > 750 and max_wind_pres <=800:
##                                height_index = 'orangered'
##                            elif max_wind_pres > 800 and max_wind_pres <=850:
##                                height_index = 'mediumvioletred'
##                            elif max_wind_pres > 850 and max_wind_pres <=900:
##                                height_index = 'darkviolet'
##                            elif max_wind_pres > 900:
##                                height_index = 'mediumblue'
#
#
### Set the map bounds
##ax.set_xlim([-66,-60])
##ax.set_ylim([-34,-28])
#
#print('saving')
#
#plt.savefig('/home/disk/meso-home/crs326/Documents/Research/WRF_Paper/PNNL_WRF/WRF_MAPS/SALLJ_frequency_maps/NEW_v4PNNL_WRF_%s_%s_relaxed_SALLJ_frequency.png' %(str(start_date)[4:], str(end_date)[4:]), dpi=200)
#
#print('saved')
##print '---total time: %s seconds ---' %(time.time() - start_time)
#
