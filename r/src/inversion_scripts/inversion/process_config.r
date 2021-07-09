## This script sets the plotting options for the analysis wrapper 
## script. 
## Created by LTK and modified by DVM on 4/23/2019

###########################################################################
#############            General paths/options               ##############  
###########################################################################
#Diagnostic plot output path?
img_path = './plots/'

#Footprint path
stat_foot_dir <- "../footprints/"

#Code to source?
source("./process_output/src/make_netcdf.r")
source("./process_output/src/make_kmz.r")

flux_downscale_time_func_script <- "../include/downscale_time_func.r"
flux_downscale_space_func_script <- "../include/downscale_space_func.r"

#Subset flux times for afternoon bins only ("flux_times"). These bins are 
#time-centered, i.e, 21z, 03z, 09z, and 15z of each day. 
iaft = seq(1, ntimes, by = 4)
###########################################################################
###########               TRAX plot options                     ###########  
###########################################################################
#What is the bin length of TRAX data points?
bin.length = 50 
