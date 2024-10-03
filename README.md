# WRF_Upscale_Growth_Paper

Repository relating to the final chapter of my Ph.D. work using the PNNL WRF run to understand environments where storms grow upscale (organize into larger systems) in Central Argentina. This chapter expands upon and confirms conclusions from case studies of observed upscale growth.

## Data
Note: data is not in respository

### WRF hourly 3D Output 
Variables: 
- height
- pressure
- temp
- relative humidity
- wind speed
- wind direction
- u wind
- v wind

### WRF MCS Tracks Output
Stats files: mcs_tracks_pf_20181015_20190430.nc, robust_mcs_tracks_20181015_20190430.nc

Variables:
- datetimestring (MCS tracked times, 0 time index for MCS initiation)
- meanlon (center latitude, 0 time index for MCS initiation)
- meanlat (center longitude, 0 time index for MCS initiation)
- starttrackresult (MCS start type, 10 is new tracked MCS)
- length (MCS duration)
- majoraxislength (time index 1 minus time index 0 to get axis growth)
- css_area (Area of cold cloud shield)


The files below were produced from script in Research/WRF_Paper/PNNL_WRF/WRF_allSALLJs and were copied over to this directory for use in Composite_shear_maps_by_SALLJ_presence:

3relaxed_noSALLJ_times_VMRS_20181101_20190430.dat
3relaxed_SALLJ_times_VMRS_20181101_20190430.dat
3relaxed_SALLJ_level_max_VMRS_20181101_20190430.dat
