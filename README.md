# WRF_Upscale_Growth_Paper

Repository relating to the final chapter of my Ph.D. work using the PNNL WRF run to understand environments where storms grow upscale (organize into larger systems) in Central Argentina. This chapter expands upon and confirms conclusions from case studies of observed upscale growth.

## Data
Not in respository

WRF hourly 3D variables 
- height, pressure, wind speed, wind direction
WRF MCS tracks (stats files: mcs_tracks_pf_20181015_20190430.nc)
- datetimestring (at 0 time for MCS initiation datetime)
- pf_lon (center latitude)
- pf_lat (center longitude)


The files below were produced from script in Research/WRF_Paper/PNNL_WRF/WRF_allSALLJs and were copied over to this path:

3relaxed_noSALLJ_times_VMRS_20181101_20190430.dat
3relaxed_SALLJ_times_VMRS_20181101_20190430.dat
3relaxed_SALLJ_level_max_VMRS_20181101_20190430.dat
