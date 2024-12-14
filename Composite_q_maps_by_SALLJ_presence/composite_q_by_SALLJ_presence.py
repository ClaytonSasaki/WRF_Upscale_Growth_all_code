#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on wed Oct 9 2024

@author: crs326

This script generates .dat files for plotting composite specific humidity (q) for SALLJ or non-SALLJ times. 
Note that this script does not produce plots; for plotting, refer to a separate script.

It is based on 'composite_map_shear_SALLJ_times.py' in Composite_shear_maps_by_SALLJ_presence directory.

Originally based on 'relaxed_SALLJ_spatial_map_frequency_efficient_v4.py', 
'plot_only_SALLJ_spatial_map_frequency_efficient_v4.py', 'plot_only_relaxed_SALLJ_distributions_by_height.py', 
and other related files.

It uses the '3relaxed_SALLJ_times_VMRS_1101_0430.dat' file. Please refer to the README for more information.

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
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap

arr = np.linspace(0, 50, 100).reshape((10, 10))
fig, ax = plt.subplots(ncols=2)

cmap = plt.get_cmap('terrain')
new_cmap = truncate_colormap(cmap, 0.21, 1)
ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

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

data_path = '/home/disk/meso-home/crs326/Documents/Research/WRF_Upscale_Growth_Paper/Composite_q_maps_by_SALLJ_presence/data/'
path_wrf = '/home/disk/monsoon/relampago/raw/wrf/'

# change the start and end days to be plotted
start_date = 20181101
end_date = 20190430

# convert dates to datetime objects to only run for correct day folders
input_start_day = pd.to_datetime(str(start_date), format='%Y%m%d', errors='ignore')
input_end_day = pd.to_datetime(str(end_date), format='%Y%m%d', errors='ignore')

# define SALLJ presence status
SALLJ_presence = 'noSALLJ' # SALLJ or noSALLJ
SALLJ_noSALLJ_location = 'VMRS'

# Define the pressure levels for analysis
q_level1 = 500
q_level2 = 700
q_level3 = 850

# load the time list for relaxed SALLJ criteria
time_list = pickle.load(open(data_path + "3relaxed_%s_times_%s_%s_%s.dat" %(SALLJ_presence, SALLJ_noSALLJ_location, start_date, end_date), "rb")) # This files used relaxed SALLJ criteria that allows for jets centered higher in the atmosphere, and a bit weaker compared to mean flow. max wind up to 3200m, minimum up to 5700m, max wind > 19.4384, decrease from max to min > 7.77538 (relaxed to 10 m/s max and 4 m/s decrease)

time_list_arr = np.array([str(time) for time in time_list])

time_list_arr_dt = [datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S') for date in time_list_arr]

# initialize empty arrays to store q values at different pressure levels
q_level1_spatial_times = np.empty((len(time_list_arr_dt), 499, 599))
q_level1_spatial_times[:] = np.NaN

q_level2_spatial_times = np.empty((len(time_list_arr_dt), 499, 599))
q_level2_spatial_times[:] = np.NaN

q_level3_spatial_times = np.empty((len(time_list_arr_dt), 499, 599))
q_level3_spatial_times[:] = np.NaN

count = 0

# ----------------- Process NetCDF Files -------------------------

for time in time_list_arr_dt:

    dt_date_str = time.date().strftime('%Y%m%d')

    file_name = '/wrfout_d01_%s-%s-%s_%s:00:00' %(time.strftime('%Y'), time.strftime('%m'), time.strftime('%d'), time.strftime('%H'))

    file_path = path_wrf + dt_date_str + file_name

    print(file_path)

    # get the netCDF
    ncfile = Dataset(file_path,'r')

    ############################# read in some variables and calculate at chosen levels #############################
    print('reading in variables for: %s %s'  %(dt_date_str, time.strftime('%H')))
    # read in the variables 
    pres = getvar(ncfile, 'pressure')
    temp_c = getvar(ncfile, 'temp', units='degC')
    RH = getvar(ncfile, 'rh') # in kg/kg

    #### interpolate to the set level1 and calculate q ####
    temp_c_level1 = np.squeeze(interplevel(temp_c, pres, q_level1))

    RH_level1 = np.squeeze(interplevel(RH, pres, q_level1))

    e_sat_level1 = SatVap(temp_c_level1)

    e_level1 = (RH_level1/100.) * e_sat_level1

    w_level1 = MixRatio(e_level1,q_level1*100.)
    #w = (e*Rs_da) / (Rs_v*(pres*100. - e))
    q_level1 = (w_level1 / (w_level1+1))*1000 #g/kg

    #### interpolate to the set level2 and calculate q ####
    temp_c_level2 = np.squeeze(interplevel(temp_c, pres, q_level2))

    RH_level2 = np.squeeze(interplevel(RH, pres, q_level2))

    e_sat_level2 = SatVap(temp_c_level2)

    e_level2 = (RH_level2/100.) * e_sat_level2

    w_level2 = MixRatio(e_level2,q_level2*100.)
    #w = (e*Rs_da) / (Rs_v*(pres*100. - e))
    q_level2 = (w_level2 / (w_level2+1))*1000 #g/kg

    #### interpolate to the set leve21 and calculate q ####
    temp_c_level3 = np.squeeze(interplevel(temp_c, pres, q_level3))

    RH_level3 = np.squeeze(interplevel(RH, pres, q_level3))

    e_sat_level3 = SatVap(temp_c_level3)

    e_level3 = (RH_level3/100.) * e_sat_level3

    w_level3 = MixRatio(e_level3,q_level3*100.)
    #w = (e*Rs_da) / (Rs_v*(pres*100. - e))
    q_level3 = (w_level3 / (w_level3+1))*1000 #g/kg

    # store the calculated q values for each time step
    q_level1_spatial_times[count,:,:] = q_level1
    q_level2_spatial_times[count,:,:] = q_level2
    q_level3_spatial_times[count,:,:] = q_level3

    count = count + 1
    ncfile.close()

# calculate mean and median over all times and save as .dat files
mean_q_level1_spatial_times = np.nanmean(q_level1_spatial_times, axis=0)
print(mean_q_level1_spatial_times)
median_q_level1_spatial_times = np.nanmedian(q_level1_spatial_times, axis=0)
print(median_q_level1_spatial_times)

pickle.dump(mean_q_level1_spatial_times, open(data_path + "mean_%s_%s_times_%s_%s.dat" %('500hPa_q', SALLJ_presence, start_date, end_date), "wb"))
pickle.dump(median_q_level1_spatial_times, open(data_path + "median_%s_%s_times_%s_%s.dat" %('500hPa_q', SALLJ_presence, start_date, end_date), "wb"))

mean_q_level2_spatial_times = np.nanmean(q_level2_spatial_times, axis=0)
print(mean_q_level2_spatial_times)
median_q_level2_spatial_times = np.nanmedian(q_level2_spatial_times, axis=0)
print(mean_q_level2_spatial_times)

pickle.dump(mean_q_level2_spatial_times, open(data_path + "mean_%s_%s_times_%s_%s.dat" %('700hPa_q', SALLJ_presence, start_date, end_date), "wb"))
pickle.dump(median_q_level2_spatial_times, open(data_path + "median_%s_%s_times_%s_%s.dat" %('700hPa_q', SALLJ_presence, start_date, end_date), "wb"))

mean_q_level3_spatial_times = np.nanmean(q_level3_spatial_times, axis=0)
print(mean_q_level3_spatial_times)
median_q_level3_spatial_times = np.nanmedian(q_level3_spatial_times, axis=0)
print(mean_q_level3_spatial_times)

pickle.dump(mean_q_level3_spatial_times, open(data_path + "mean_%s_%s_times_%s_%s.dat" %('850hPa_q', SALLJ_presence, start_date, end_date), "wb"))
pickle.dump(median_q_level3_spatial_times, open(data_path + "median_%s_%s_times_%s_%s.dat" %('850hPa_q', SALLJ_presence, start_date, end_date), "wb"))