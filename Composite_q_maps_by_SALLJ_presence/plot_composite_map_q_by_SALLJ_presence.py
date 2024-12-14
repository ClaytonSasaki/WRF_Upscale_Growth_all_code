#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Oct 9 2024

@author: crs326

Plots composite maps of q for all SALLJ times or NON-SALLJ times. Takes data created by 'composite_q_by_SALLJ_presence.py' script

"""

import matplotlib
matplotlib.use('Agg') # For non-interactive plotting (e.g., on a server)

import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature, LAND, OCEAN, COASTLINE, BORDERS, LAKES, RIVERS
import datetime
import pickle
import math
import matplotlib.colors as colors
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
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, ll_to_xy, get_basemap, xy_to_ll, CoordPair, GeoBounds)

# ------------------------- Custom Colormaps ------------------------------

# used for contourf of terrain 
def truncate_colormap2(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap2 = plt.get_cmap('gray')
new_cmap2 = truncate_colormap(cmap2, 0.3, 0.75)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap2)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap2)

# ----------------- Define Atmospheric Calculations ------------------

def SatVap(tempc,phase="liquid"):
    """Calculate saturation vapour pressure over liquid water and/or ice.

    INPUTS: 
    tempc: (C)
    phase: ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
    return saturation vapour pressure as follows:

    Tc>=0: es = es_liquid
    Tc <0: es = es_ice


    RETURNS: e_sat  (Pa)

    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)

    This formulation is chosen because of its appealing simplicity, 
    but it performs very well with respect to the reference forms
    at temperatures above -40 C. At some point I'll implement Goff-Gratch
    (from the same resource).
    """

    over_liquid=6.112*exp(17.67*tempc/(tempc+243.12))*100.
    over_ice=6.112*exp(22.46*tempc/(tempc+272.62))*100.
    # return where(tempc<0,over_ice,over_liquid)

    if phase=="liquid":
    # return 6.112*exp(17.67*tempc/(tempc+243.12))*100.
        return over_liquid
    elif phase=="ice":
    # return 6.112*exp(22.46*tempc/(tempc+272.62))*100.
        return where(tempc<0,over_ice,over_liquid)
    else:
        raise NotImplementedError

Epsilon=0.622         # Epsilon=R_s_da/R_s_v; The ratio of the gas constants

def MixRatio(e,p):
    """Mixing ratio of water vapour
    INPUTS
    e (Pa) Water vapor pressure
    p (Pa) Ambient pressure

    RETURNS
    qv (kg kg^-1) Water vapor mixing ratio`
    """

    return Epsilon*e/(p-e)

Rs_da=287.05          # Specific gas const for dry air, J kg^{-1} K^{-1}
Rs_v=461.51           # Specific gas const for water vapour, J kg^{-1} K^{-1}

# ----------------- File Paths and Settings -------------------------

general_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/Composite_q_maps_by_SALLJ_presence/'
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

# change the start and end days to be plotted
start_date = 20181101
end_date = 20190430

# convert dates to datetime objects to only run for correct day folders
input_start_day = pd.to_datetime(str(start_date), format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(end_date), format='%Y%m%d', errors='ignore')

# define SALLJ presence status
SALLJ_presence = 'SALLJ' # SALLJ or noSALLJ
SALLJ_noSALLJ_location = 'VMRS'

# define plotted variables
plotted_var = '850hPa_q'
mean_median = 'mean'

################################# plot #################################
terrain = pickle.load(open(general_path + "data/terrain.dat", "rb"))
lats = pickle.load(open(general_path + "data/lats.dat", "rb"))
lons = pickle.load(open(general_path + "data/lons.dat", "rb"))
cart_proj = pickle.load(open(general_path + "data/cart_proj.dat", "rb"))

# Create the figure
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

# Plot markers for COR and VMRS
ax.scatter(-64.212, -31.298, s=30, color='black', transform=crs.PlateCarree(), zorder=6)
ax.scatter(-63.726, -29.906, s=30, color='black', transform=crs.PlateCarree(), zorder=6)

print('plotted terrain')

# Load file with variable data to plot
mean_median_variable_spatial_SALLJ_times = pickle.load(open(general_path + "data/%s_%s_SALLJ_times_%s_%s.dat" %(mean_median, plotted_var, start_date, end_date), "rb"))

# Smooth variable data
mean_median_variable_spatial_SALLJ_times_smooth = smooth2d(mean_median_variable_spatial_SALLJ_times, 6, cenweight=4)

# created mask where terrain is less than 3000m and only use data from that area
mask = terrain <= 3000
mean_median_variable_spatial_SALLJ_times_smooth_terrainLess3000 = np.where(mask, mean_median_variable_spatial_SALLJ_times_smooth, np.nan)

# plot variable data contoured
levels = np.arange(0,20,2)

var_plot = plt.contour(lons, lats, mean_median_variable_spatial_SALLJ_times_smooth_terrainLess3000, levels=levels, cmap='Reds', transform=crs.PlateCarree())

ax.clabel(var_plot, var_plot.levels, inline=True, fmt='%.0f')

plt.tight_layout()

print('saving')

plt.savefig(general_path + 'plots/%s_%s_composite_map_%s_times_%s_%s_%s.png' %(plotted_var, mean_median, SALLJ_presence, start_date, end_date, SALLJ_noSALLJ_location))

print('saved')