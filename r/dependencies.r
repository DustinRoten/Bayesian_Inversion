### Source R scripts for Dien's X-STILT, 05/23/2018, DW
rsc <- dir('r/src', pattern = '.*\\.r$',
           full.names = T, recursive = F)

invisible(lapply(rsc, source))

# Other relevant packages
library(ggplot2); library(viridis); library(patchwork)
library(geosphere); library(geodist); library(raster)
library(ggmap); library(ncdf4); library(dplyr)
library(lutz); library(lubridate); library(stringr)
library(reshape2)

# required STILT functions
#source('X-STILT/stilt/r/src/write_footprint.r')

# required X-STILT functions
#source('X-STILT/r/src/get.lon.lat.r')
