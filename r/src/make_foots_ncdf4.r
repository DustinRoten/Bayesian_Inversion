make_foots_ncdf4 <- function(footprint.list = NULL, averaged.timesteps = NULL,
                             timestr = NULL, operator.dir = NULL, nc.extent = NULL) {
  
  for(i in 1:length(footprint.list)) {
    
    #output message
    cat(paste0('Footprint formatting ',
               round(100*i/length(footprint.list), 2),
               '% complete.     '), '\r')
    
    ##########################################
    ##### Read in and crop the footprint #####
    ##########################################
    #create the footprint's personal directory
    SAM.name <- paste0('footprint_', i)
    if(!dir.exists(file.path(operator.dir, SAM.name)))
      dir.create(file.path(operator.dir, SAM.name),
                 recursive = TRUE)
    
    #get filename and create the full file name for new footprint
    file.name <- gsub('_X_foot.nc', '',
                      basename(footprint.list[i]))
    
    #replace the old date/time from X-STILT with the shifted date/time
    #first, identify the components making up the file name
    file.name.components <- unlist(strsplit(file.name, split = '_'))
    
    #which one matches?
    partial.match <- grepl(pattern = substr(timestr, 1, 8),
                           file.name.components)
    part.idx <- which(partial.match)
    
    #construct the new name using the averaged start time
    file.name.components[part.idx] <- timestr
    new.file.name <- paste(file.name.components, collapse = '_')
      
    #finally, generate the full name and path of the file
    file.fullname <- file.path(operator.dir, SAM.name,
                               paste0(new.file.name, '.nc'))
    
    #read in and crop the footprint
    footprint <- brick(footprint.list[i])
    cropped.footprint <- crop(footprint, nc.extent)
    remove('footprint') #no longer need this (save RAM!)
    
    #' ncdf metadata is NOT saved in `writeRaster()`
    layer.names <- gsub('X', '', names(cropped.footprint))
    layer.names.pos <- as.POSIXlt(layer.names, tz = 'UTC',
                                  format = '%Y.%m.%d.%H.%M.%S')
    
    ##################################
    ##### Update the layer names #####
    ##################################
    new.names <- NULL
    for(j in 1:nlayers(cropped.footprint)) {
      
      #map the layer timestamp to one of the averaged timesteps
      #(time difference in seconds)
      name.idx <-
        which(abs(as.numeric(averaged.timesteps) -
                    as.numeric(layer.names.pos[j])) < 380)
      
      #' If there are no matches, then there is likely an issue
      #' with the calculated averaged timesteps.
      if(length(name.idx) == 0 | length(name.idx) > 1)
        stop(paste0('Error in timesteps: ', '\n',
                    'Footprint: ', footprint.list[i], '\n',
                    'Iteration: ', j, '\n',
                    'Error Type: ', length(name.idx)))
      
      #add the averaged timestep to the vector
      new.names[j] <- as.numeric(averaged.timesteps[name.idx])
      
    }; names(cropped.footprint) <- new.names
    
    # save the raster here, reopen and update information
    suppressWarnings(writeRaster(cropped.footprint,
                                 filename = file.fullname,
                                 format = 'CDF', overwrite = TRUE))
    remove('cropped.footprint')
  
    #read the file in again as a NetCDF file
    nc_f <- nc_open(file.fullname, write = TRUE)
    
    #' Define these dimensions as prescribed by L. Kunik and D. Mallia
    #' in previous code.
    #convert time
    time_dim <-
      ncdim_def("time", "seconds_since_1970_01_01",
                new.names,
                longname = "seconds since R epoch: 1970-01-01 00:00:00")
    ntime <- length(time_dim$vals)
    
    #convert lats
    lat_dim <- ncdim_def("lat", "degrees_north",
                         ncvar_get(nc_f, 'latitude'),
                         longname = "latitude (center of cell)")
    nlat <- length(lat_dim$vals)
    
    #convert lons
    lon_dim <- ncdim_def("lon", "degrees_east",
                         ncvar_get(nc_f, 'longitude'),
                         longname = "longitude (center of cell)")
    nlon <- length(lon_dim$vals)
    
    #prepare the footprint variable
    foot_var <-
      ncvar_def("foot", "ppm/(umol m-2 s-1)",
                list(lon_dim, lat_dim, time_dim),
                longname = "gridded footprint influence values (mole fraction per flux)")
    
    foot_vars <- list(foot_var)
    
    foot_values <- ncvar_get(nc_f, "time")
    
    #' default name after `writeRaster()` is 'time'. This must
    #' be changed to the proper footprint value (Kunik & Mallia).
    nc_filename <- file.fullname
    nc_newfoot <- nc_create(nc_filename, foot_vars)
    ncvar_put(nc_newfoot, foot_var, foot_values, count = c(nlon, nlat, ntime))
    nc_close(nc_newfoot)
    
  }; cat('') #new line
  
}