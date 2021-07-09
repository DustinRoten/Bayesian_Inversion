#Compile biogenic fluxes into inversion input format. If we are reading into more than one 
#monthly file at a time, concatenate these netcdf files. Note that this code is currently 
#written to work with Doug Catherine's biospheric flux code, and will need to be adapted to 
#work for other biospheric flux inventories.
#
# Written by Lewis Kunik, modified by DVM 6/3/2019
# Adapted for SMUrF by Dustin Roten 5/19/2021

#' In cases where you require data from an incomplete SMUrF
#' dataset, include an alternative year as `match.year` to
#' pull biospheric flux from.
#' 
#' Additionally, supply any raster to `match.raster`. This is the
#' raster that SMUrF will be cropped and resampled onto.

make_bio_ncdf4 <- function(SMUrF.path = NULL, times = NULL,
                           match.year = NA, match.raster = NULL,
                           bio.output = NULL) {
  
  #require appropriate libraries
  require(ncdf4)
  require(raster)
  require(abind)
  
  #make sure the times vector is properly formatted
  if(is.POSIXt(times) == FALSE) stop('Times are not in POSIX format!')
  
  #############################################################
  #' If `match.year == T`, then a previous year will be used if
  #' current year is not available.
  #############################################################
  all.SMUrF.files <- list.files(SMUrF.path, full.names = TRUE)
  if(is.na(match.year)) {
  
    #change the times to YYYYMM format to match SMUrF output files
    year.month <- unique(strftime(times, format = '%Y%m'))
    
    #Get the indices of corresponding SMUrF files.
    matched.files <-
      which(!is.na(str_match(basename(all.SMUrF.files),
                             pattern = year.month)))
    SMUrF.files <- all.SMUrF.files[matched.files]
    
    #If there aren't enough SMUrF files, stop.
    if(length(year.month) > length(SMUrF.files))
      stop('Missing SMUrF files! Consider using match.year == F.')
    
  } else if(!is.na(match.year)) {
    
    message('Using alternative year for biospheric flux.')
    
    #change the times to YYYYMM format to match SMUrF output files
    year.month <- unique(strftime(times, format = '%Y%m'))
    
    #replace the current year with the desired `match.year`
    year.month <- gsub(unique(strftime(times, format = '%Y')),
                       paste0(match.year), year.month)
    
    #Get the indices of corresponding SMUrF files.
    matched.files <-
      which(!is.na(str_match(basename(all.SMUrF.files),
                             pattern = year.month)))
    SMUrF.files <- all.SMUrF.files[matched.files]
  }
  
  #' Construct a lookup table for each SMUrF file and its 
  #' associated hourly timesteps.
  for(i in 1:length(SMUrF.files)) {
    
    #read in the SMUrF raster
    tmp.raster <- brick(SMUrF.files[i])
    layer.names <- names(tmp.raster)
    
    #add to the lookup table
    if(i == 1)
      lookup.table <- 
      data.frame(file = SMUrF.files[i],
                 layer.names = layer.names,
                 POSIX = as.POSIXlt(gsub('X', '', layer.names),
                                    format = '%Y.%m.%d.%H.%M.%S',
                                    tz = 'UTC'))
    if(i > 1)
      lookup.table <-
      rbind(lookup.table,
            data.frame(file = SMUrF.files[i],
                       layer.names = layer.names,
                       POSIX = as.POSIXlt(gsub('X', '', layer.names),
                                          format = '%Y.%m.%d.%H.%M.%S',
                                          tz = 'UTC')))
    
  }; remove('tmp.raster')
  ############################################################
  
  
  #################################################
  #' Begin looping through the `times` and making a
  #' raster stack of SMUrF output.
  #################################################
  
  for(i in 1:length(times)) {
    
    #' `tmp.time` and `lookup.table` may be off by a year.
    #' This is remedied here by restructuring `tmp.time`.
    if(!is.na(match.year)) {
      
      #create a new time to the nearest hourly resolution
      #strip the minutes and hours off
      tmp.time <- as.numeric(times[i])
      tmp.time <- tmp.time - 60*minute(times[i])
      new.time <- tmp.time - second(times[i])
      
      #also have new.time as a pos.time
      tmp.pos.time <- 
        as.POSIXlt(new.time, origin = '1970-01-01', tz = 'UTC')
      
      #create a temporary format to be read in as a POSIX time
      match.time <- paste(match.year,
                          month(tmp.pos.time),
                          day(tmp.pos.time),
                          hour(tmp.pos.time),
                          sep = '-')
      match.time <- as.POSIXlt(match.time,
                               format = '%Y-%m-%d-%H',
                               tz = 'UTC')
      
      #find this time in the lookup table
      tbl.idx <- which(lookup.table$POSIX == match.time) 
      
    } else if(is.na(match.year)) {
      
      #if match.year == NA, then there should be an exact match
      which(lookup.table$POSIX == match.time)
      tbl.idx <- which(lookup.table$POSIX == match.time)
      
    }
    
    #' read in the raster layer, crop it, and modify the resolution
    #' to match the outer emissions domain.
    smurf.raster <- brick(lookup.table$file[tbl.idx])
    smurf.raster.layer <-
      eval(parse(text = paste0('smurf.raster$',
                               lookup.table$layer.names[tbl.idx])))
    remove('smurf.raster') #no longer needed
    
    #resample and crop the bio layer. All NA's are set to zero.
    smurf.raster.layer <-
      resample(smurf.raster.layer, match.raster, method = 'ngb')
    smurf.raster.layer <- crop(smurf.raster.layer, match.raster)
    smurf.raster.layer[is.na(smurf.raster.layer)] <- 0
    
    if(i == 1) bio.layers <- smurf.raster.layer
    if(i > 1) bio.layers <- addLayer(bio.layers, smurf.raster.layer)
    
  }; remove('smurf.raster.layer')
  
  #once all of the layers have been added, save the stack
  bio.layers[is.na(bio.layers)] <- 0
  
  #correct times here if an alternative year was used.
  if(!is.na(match.year)) {
    
    #' Grab the `times` and convert them to epoch time.
    #' This new time format will be used as the layer
    #' names in the biospheric flux file. (SMURF)
    names(bio.layers) <- paste0('X', as.numeric(times))
     
  } else if(is.na(match.year)) {
    
    #read in the character string, convert to POSIX, paste as numeric
    converted.layer.names <-
      as.POSIXlt(gsub('X', '', names(bio.layers)),
                 format = '%Y.%m.%d.%H.%M.%S',
                 tz = 'UTC')
    names(bio.layers) <- as.numeric(converted.layer.names)
    
  }
  
  #save the raster
  writeRaster(bio.layers, filename = bio.output,
              format = 'CDF', overwrite = TRUE)
  remove('bio.layers')
  
  #OPEN AS A NetCDF FILE AND CORRECT THE DIM/VAR FORMATS
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(bio.output, write = TRUE)
  
  #' Define these dimensions as prescribed by L. Kunik and D. Mallia
  #' in previous code.
  #convert time
  time_dim <- ncdim_def("time", "seconds_since_1970_01_01",
                        as.numeric(times),
                        longname = "seconds since R epoch: 1970-01-01 00:00:00")
  ntime <- length(time_dim$vals)
  
  #convert lats
  lat_dim <- ncdim_def("lat", "degrees_north",
                       ncvar_get(nc, 'latitude'),
                       longname = "latitude (center of cell)")
  nlat <- length(lat_dim$vals)
  
  #convert lons
  lon_dim <- ncdim_def("lon", "degrees_north",
                       ncvar_get(nc, 'longitude'),
                       longname = "longitude (center of cell)")
  nlon <- length(lon_dim$vals)  

  #assign these new dimensions and variables
  bio_flux_var <- ncvar_def("flux", "umol m-2 s-1",
                            list(lon_dim, lat_dim, time_dim),
                            longname = "bio fluxes")
  bio_vars <- list(bio_flux_var)
  
  #get the variable values from the raster.
  #(defaults to "variable")
  bio_values <- ncvar_get(nc, "variable")
  nc_close(nc)
  
  #' Here, default name after `writeRaster()` is 'layer'.
  #' This must be changed to the proper bio flux value 
  #' (Kunik & Mallia).
  nc_bio <- nc_create(bio.output, bio_vars)
  ncvar_put(nc_bio, bio_flux_var, bio_values, count = c(nlon, nlat, ntime))
  nc_close(nc_bio)

}