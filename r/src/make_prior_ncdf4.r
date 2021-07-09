make_prior_ncdf4 <- function(site = NULL, odiac.path = NULL, edgar.path = NULL,
                             category = NA, sector.list = NULL, times = NULL,
                             inner.extent = NULL, full.extent = NULL,
                             downscaling.extension = NULL, carma.path = NULL,
                             obs.dir = NULL, inner.name = NULL, outer.name = NULL) {
  
  #check the format of the input times
  if(!any(class(times) == 'POSIXt'))
    stop('Input times must be POSIXt format!')
  
  require(ggmap)
  lon.lat <- 
    suppressWarnings(get.lon.lat(site,
                                 dlon = (inner.extent[2]-inner.extent[1])/2,
                                 dlat = (inner.extent[4]-inner.extent[3])/2))
                    
  prev.YYYYMM <- 'NA'
  for(i in 1:length(times)) {
    
    #format the timestep
    YYYYMM <- strftime(times[i], format = '%Y%m')
    
    #get ODIAC layer
    if(YYYYMM != prev.YYYYMM) {
      ODIAC <- get.odiac(tiff.path = odiac.path,
                         nc.extent = full.extent,
                         YYYYMM, convert.units = FALSE)
      
      #identify large point sources
      ODIAC.lps <-
        get.odiac.lps(CarMA = read.csv(carma.path), ODIAC)
      
      #identify areal sources
      ODIAC.nlps <-
        get.odiac.nlps(CarMA = read.csv(carma.path), ODIAC)
      
    }; prev.YYYYMM <- YYYYMM
    
    #determine total.edgar for the timestep here
    total.edgar <- 
      timestep.normalization.constant(lon.lat,
                                      sector_df = sector.list,
                                      ODIAC,
                                      carma.inventory = carma.path,
                                      edgar.inventory = edgar.path,
                                      nc.extent = full.extent,
                                      temporal.downscaling = 
                                        downscaling.extension,
                                      time =
                                        strftime(times[i],
                                                 format = '%Y.%m.%d.%H.%M.%S',
                                                 tz = 'UTC'))
    
    #for now, only the inner domain will be temporally downscaled.
    #get all of the sectors
    sectors <- NULL
    if(is.na(category)) {
      
      for(j in 1:nrow(sector.list)) {
        tmp <-
          unlist(strsplit(sector.list$EDGAR_SectorList[j], split = ' '))
        sectors <- c(sectors, tmp)
      }
      
    } else if(!is.na(category)) {
      tmp <-
        unlist(strsplit(sector.list$EDGAR_SectorList[sector.list == category],
                        split = ' '))
      sectors <- c(sectors, tmp)
    }
    ##############################
    ##############################
    ##############################
    
    
    ######################################
    ##### Start weighting ODIAC here #####
    ######################################
    for(j in 1:length(sectors)) {
      
      if(sectors[j] != 'ENE') {
        
        #Grab the edgar sector (that is not ENE)
        #Resample this sector at ODIAC resolution
        edgar.sector <-
          get.edgar.sector(nc.path = list.files(edgar.path,
                                                pattern = sectors[j],
                                                full.names = TRUE),
                           nc.extent = full.extent)
        edgar.sector <- resample(edgar.sector, ODIAC)
        edgar.sector[is.na(edgar.sector)] <- 0
        edgar.sector[edgar.sector < 0] <- 0
        
        #get the appropriate weighting for the sector
        weight <- 
          edgar.sector.weighting(citylon = lon.lat$citylon,
                                 citylat = lon.lat$citylat,
                                 local.tz = lon.lat$tz,
                                 sector.name = sectors[j],
                                 temporal.downscaling.files = 
                                   downscaling.extension,
                                 time =
                                   strftime(times[i],
                                            format = '%Y.%m.%d.%H.%M.%S',
                                            tz = 'UTC'),
                                 monthly = TRUE)
        
        #weight the edgar sector/layer
        # This EDGAR raster is in kg / m2 / hr.
        # convert the unit of CO2 emiss from kg/m2/hr to umol/m2/s
        edgar.sector <- edgar.sector * (10^9 / 44.01) # convert kg to uomol-C (= umole-CO2)
        edgar.sector <- edgar.sector / 60 / 60 # convert per hour to per second
        weighted.ODIAC <- (weight*edgar.sector/total.edgar)*ODIAC.nlps
        # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
        
        #' `weighted.ODIAC` still has the original ODIAC units: Tonnes/cell/month.
        #' These need to be converted and saved.
        #compute area using area() function in raster package
        area.raster <- raster::area(weighted.ODIAC) * 1E6    # convert km2 to m2
        
        # how many days are in the particular month?
        month.days <- days_in_month(times[i])
        
        # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
        weighted.ODIAC <- weighted.ODIAC * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
        weighted.ODIAC <- weighted.ODIAC / month.days / 24 / 60 / 60	# convert per month to per second
        weighted.ODIAC <- weighted.ODIAC / area.raster		# convert per cell to per m2
        # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
        remove('area.raster') #save RAM
        
        #add the edgar layer here
        if(j == 1)
          edgar.layers <- weighted.ODIAC
        if(j > 1)
          edgar.layers <- edgar.layers + weighted.ODIAC
        
        remove('weighted.ODIAC') #save RAM
        
      } else if(sectors[j] == 'ENE') { #closes !ENE
        
        #get the sector weighting
        weight <-
          edgar.sector.weighting(citylon = lon.lat$citylon,
                                 citylat = lon.lat$citylat,
                                 local.tz = lon.lat$tz,
                                 sector.name = sectors[j],
                                 temporal.downscaling.files =
                                   downscaling.extension,
                                 time = strftime(times[i],
                                                 format = '%Y.%m.%d.%H.%M.%S',
                                                 tz = 'UTC'),
                                 monthly = TRUE)
        
        #construct with the LPS version of ODIAC
        pwr.sector <- weight*ODIAC.lps
        
        # multiplying by the temporal weight reduces it to Tonnes / cell / hr
        # compute area using area() function in raster package
        area.raster <- raster::area(pwr.sector) * 1E6    # convert km2 to m2
        
        # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
        pwr.sector <- pwr.sector * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
        pwr.sector <- pwr.sector / 60 / 60	# convert per hour to per second
        pwr.sector <- pwr.sector / area.raster		# convert per cell to per m2
        # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
        remove('area.raster') #save RAM
        
        #add to the layers
        if(j == 1) edgar.layers <- pwr.sector
        if(j > 1) edgar.layers <- edgar.layers + pwr.sector
        remove('pwr.sector') #save RAM
        
      } #close ENE
      
    } #close sector loop
    
    #add inner and outer layers to the output stacks
    if(i == 1) {
      #save the inner extent
      downscaled.ODIAC.inner <- crop(edgar.layers, inner.extent)
      
      #save the outer extent
      #create the null inner area
      downscaled.ODIAC.outer <-
        edgar.layers*extend(0*downscaled.ODIAC.inner,
                            full.extent, value = 1)
    } else if(i > 1) {
      #if the process has already started...
      downscaled.ODIAC.inner <-
        addLayer(downscaled.ODIAC.inner,
                 crop(edgar.layers, inner.extent))
      
      #save the outer extent
      #create the null inner area
      #' This command will throw a warning due to the mismatched
      #' domains of the inner, outer, and EDGAR rasters. This is
      #' fine... cropping will happen automatically in the
      #' `AddLayer` function. Discrepancies can be handled later.
      downscaled.ODIAC.outer <- suppressWarnings(
        addLayer(downscaled.ODIAC.outer,
                 edgar.layers*extend(0*downscaled.ODIAC.inner[[1]],
                                     full.extent, value = 1))
      )
    }
    
    #include a progress message
    message <- paste0('\n', 'Calculated prior emissions for time ',
                      times[i], '. ', round(100*i/length(times), 2),
                      '% complete.')
    cat(message, '\r')
    
  } #closes time loop
  
  ########################################
  ##### Save the prior_emiss.nc File #####
  ########################################
  #once all of the layers have been added, save the stack
  downscaled.ODIAC.inner[is.na(downscaled.ODIAC.inner)] <- 0
  writeRaster(downscaled.ODIAC.inner,
              filename = file.path(obs.dir, inner.name),
              format = 'CDF', overwrite = TRUE)
  
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(file.path(obs.dir, inner.name), write = TRUE)
  
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
  prior_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                               list(lon_dim, lat_dim, time_dim),
                               longname = "prior emissions")
  prior_vars <- list(prior_emiss_var)
  
  #get the variable values from the raster.
  #(defaults to "variable")
  prior_values <- ncvar_get(nc, "variable")
  nc_close(nc)
  
  #' Here, default name after `writeRaster()` is 'layer'.
  #' This must be changed to the proper footprint value 
  #' (Kunik & Mallia).
  nc_filename <- file.path(obs.dir, inner.name)
  nc_emiss <- nc_create(nc_filename, prior_vars)
  ncvar_put(nc_emiss, prior_emiss_var, prior_values, count = c(nlon, nlat, ntime))
  nc_close(nc_emiss)
  ########################################
  ########################################
  ########################################
  
  
  ########################################
  ##### Save the outer_emiss.nc File #####
  ########################################
  #once all of the layers have been added, save the stack
  downscaled.ODIAC.outer[is.na(downscaled.ODIAC.outer)] <- 0
  writeRaster(downscaled.ODIAC.outer,
              filename = file.path(obs.dir, outer.name),
              format = 'CDF', overwrite = TRUE)
  
  #open the raster as a NetCDF file and correct the information
  nc <- nc_open(file.path(obs.dir, outer.name), write = TRUE)
  
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
  outer_emiss_var <- ncvar_def("emiss", "umol m-2 s-1",
                               list(lon_dim, lat_dim, time_dim),
                               longname = "outer emissions")
  outer_vars <- list(outer_emiss_var)
  
  #get the variable values from the raster.
  #(defaults to "variable")
  outer_values <- ncvar_get(nc, "variable")
  nc_close(nc)
  
  #' Here, default name after `writeRaster()` is 'layer'.
  #' This must be changed to the proper footprint value 
  #' (Kunik & Mallia).
  nc_filename <- file.path(obs.dir, outer.name)
  nc_outer <- nc_create(nc_filename, outer_vars)
  ncvar_put(nc_outer, outer_emiss_var, outer_values, count = c(nlon, nlat, ntime))
  nc_close(nc_outer)
  ########################################
  ########################################
  ########################################
  
} #closes function