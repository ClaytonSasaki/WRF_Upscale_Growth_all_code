#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Nov 17 14:16:04 2021

@author: crs326

Based upon maps_v2.py in the folder old

Makes maps that show the frequency of when the SALLJ if found.

"""

import matplotlib
matplotlib.use('Agg') 

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import numpy as np
import numpy.ma as ma
from numpy import exp,where,ma,cos,sin,pi,amax,amin
import pandas as pd
import os
from scipy import stats
from scipy.interpolate import interp1d
import math
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import xarray
#import time
import pickle
import wrf

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

cmap = plt.get_cmap('Reds')
new_cmap = truncate_colormap(cmap, 0.0, 0.8333333333)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

# used for contour of terrain
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('gray')
new_cmap2 = truncate_colormap(cmap2, 0.4, 1.0)
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

#change the start and end days to be plotted
start_date = 20190101
end_date = 20190228

terrain = pickle.load(open("relaxed_SALLJ_frequency_v4_terrain.dat", "rb"))
lats = pickle.load(open("relaxed_SALLJ_frequency_v4_lats.dat", "rb"))
lons = pickle.load(open("relaxed_SALLJ_frequency_v4_lons.dat", "rb"))
cart_proj = pickle.load(open("relaxed_SALLJ_frequency_v4_cart_proj.dat", "rb"))

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

# set terrain plotted to True
terrain_plotted = True

print('plotted terrain')


pres_max_wind_SALLJs_for_all_times = pickle.load(open("relaxed_SALLJ_frequency_v4_%s_%s.dat" %(start_date, end_date), "rb"))

pres_max_wind_SALLJs_for_all_times_smooth = wrf.smooth2d(pres_max_wind_SALLJs_for_all_times, 6, cenweight=4)

# For 2 days (to get levels 0 to 55 by 5%)
#plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[2.4,4.8,7.2,9.6,12,14.4,16.8,19.2,21.6,24,26.4], extend='max', cmap='Greys', zorder=4, transform=crs.PlateCarree())

# For 30 days (to get levels 0 to 55 by 5%)
#plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[36,72,108,144,180,216,252,288,324,360,396], extend='max', cmap='Reds', zorder=4, transform=crs.PlateCarree())

# For 59 days (to get levels 0 to 55 by 5%)
#plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[70.8,141.6,212.4,283.2,354,424.8,495.6,566.4,637.2,708,778.8], extend='max', cmap='Reds', zorder=4, transform=crs.PlateCarree())

# For 61 days (to get levels 0 to 55 by 5%)
#plt.contourf(lons, lats, pres_max_wind_SALLJs_for_all_times, levels=[73.2,146.4,219.6,292.8,366,439.2,512.4,585.6,658.8,732,805.2], extend='max', cmap='Reds', zorder=4, transform=crs.PlateCarree())

# convert date strings to datetimes, then calculate number of hours
start_day_dt = pd.to_datetime(str(start_date), format='%Y%m%d', errors='ignore')
end_day_dt = pd.to_datetime(str(end_date), format='%Y%m%d', errors='ignore')

delta_dt = end_day_dt - start_day_dt

num_hours = (delta_dt.days+1)*24

print('num_hours', num_hours)

# choose percentage levels to plot
start_percent = 0.10 # (5 - 55 by 5 percent)
end_percent = 0.55
interval_percent = 0.05

percentages = np.arange(start_percent, end_percent, interval_percent)

levels = [i * num_hours for i in percentages]

print('percentages', percentages)

SALLJ_frequency = plt.contour(lons, lats, pres_max_wind_SALLJs_for_all_times_smooth, levels=levels, extend='max', cmap='Reds', zorder=4, transform=crs.PlateCarree())

# use for those with smaller range to keep consistent cmap
#SALLJ_frequency = plt.contour(lons, lats, pres_max_wind_SALLJs_for_all_times_smooth, levels=levels, extend='max', cmap=new_cmap, zorder=4, transform=crs.PlateCarree())

# For inline labels
#fmt = {}
#for l, s in zip(SALLJ_frequency.levels, percentages):
#    fmt[l] = '%.0f %%' %(s*100)
#
#ax.clabel(SALLJ_frequency, SALLJ_frequency.levels, inline=True, fmt=fmt)

# For legend
#for i in range(len(percentages)-1):
#    SALLJ_frequency.collections[i].set_label('%.0f %%' %(percentages[i]*100))
#
#plt.legend()

plt.tight_layout()

#ax.clabel(SALLJ_frequency, percentages, inline=True, fmt='%d %%')

## write data to .csv file
#with open('PNNL_WRF_SDC_CAPE_time_series_%s_%s_%s.csv' % (area_names[i], str(start_date)[4:], str(end_date)[4:]), 'w') as f:
#    writer = csv.writer(f)
#
#    # write the times to a row
#    writer.writerow(time_plot)
#
#    # write the data to the next rows
#    writer.writerow(Average_MUCAPE_subset_list)
#
#print('done')



#    # get the pressure of max_wind of points that meet SALLJ criteria that fall into pressure groupings
#    pres_max_wind_SALLJs_less700 = pres_max_wind_SALLJs.where(pres_max_wind_SALLJs <= 700.0, np.nan)
#    pres_max_wind_SALLJs_700_800 = pres_max_wind_SALLJs.where((pres_max_wind_SALLJs > 700.0) & (pres_max_wind_SALLJs <= 800.0), np.nan)
#    pres_max_wind_SALLJs_800_900 = pres_max_wind_SALLJs.where((pres_max_wind_SALLJs > 800.0) & (pres_max_wind_SALLJs <= 900.0), np.nan)
#    pres_max_wind_SALLJs_grt900 = pres_max_wind_SALLJs.where(pres_max_wind_SALLJs >= 900.0, np.nan)
#                            
#                            if max_wind_pres <= 700:
#                                height_index = 'gold'
#                            elif max_wind_pres > 700 and max_wind_pres <=750:
#                                height_index = 'orange'
#                            elif max_wind_pres > 750 and max_wind_pres <=800:
#                                height_index = 'orangered'
#                            elif max_wind_pres > 800 and max_wind_pres <=850:
#                                height_index = 'mediumvioletred'
#                            elif max_wind_pres > 850 and max_wind_pres <=900:
#                                height_index = 'darkviolet'
#                            elif max_wind_pres > 900:
#                                height_index = 'mediumblue'


## Set the map bounds
#ax.set_xlim([-66,-60])
#ax.set_ylim([-34,-28])

print('saving')

plt.savefig('/home/disk/meso-home/crs326/Documents/Research/WRF_Paper/PNNL_WRF/WRF_MAPS/SALLJ_frequency_maps/8v4PNNL_WRF_%s_%s_relaxed_SALLJ_frequency.png' %(str(start_date)[4:], str(end_date)[4:]), dpi=200)

print('saved')
#print '---total time: %s seconds ---' %(time.time() - start_time)

