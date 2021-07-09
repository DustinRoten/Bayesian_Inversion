synthetic.obs <- function(footprints = NULL, start.time = NULL, errors = NULL, obs.error = NA,
                          inner.emissions = NULL, outer.emissions = NULL, biosphere.emissions = NULL,
                          obs.output = NULL, background.output = NULL, std.bkg = T) {
  
  #read in the inner and outer emissions files
  inner.emissions.raster <- brick(inner.emissions)
  outer.emissions.raster <- brick(outer.emissions)
  
  #bio emissions are optional
  if(!is.na(biosphere.emissions))
    bio.emissions.raster <- brick(biosphere.emissions)
  
  site_code.list <- NULL
  XCO2_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  names(XCO2_df) <- c('date.time', 'lon', 'lat', 'xco2')
  bkg.XCO2_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 2))
  names(bkg.XCO2_df) <- c('time', 'avg_CO2')
  for(i in 1:length(footprints)) {
    
    site_code.list[i] <- paste0('footprint_', i)
    
    #begin with a progress message
    cat(paste0('Generating synthetic observations: ',
               round(100*i/length(footprints), 2),
               '% complete.     '), '\r')
    
    #get date.time here
    tmp.date.time <- 
      unlist(strsplit(basename(footprints[i]), split = '_'))[1]
    
    #get lon here
    tmp.lon <- 
      unlist(strsplit(basename(footprints[i]), split = '_'))[2]
    tmp.lon <- gsub('.nc', '', tmp.lon)
    
    #get lat here
    tmp.lat <- 
      unlist(strsplit(basename(footprints[i]), split = '_'))[3]
    tmp.lat <- gsub('.nc', '', tmp.lat)
    
    tmp.footprint <- brick(footprints[i])
    tmp.layers <- names(tmp.footprint)
    
    hourly.XCO2 <- NULL; hourly.bkg.XCO2 <- NULL
    for(j in 1:length(tmp.layers)) {
      
      tmp.footprint.timestep <- tmp.footprint[[j]]
      
      tmp.foot.layer <- gsub('X', '', tmp.layers[j])
      tmp.foot.layer.pos <- as.POSIXlt(tmp.foot.layer,
                                       format = '%Y.%m.%d.%H.%M.%S',
                                       origin = '1970-01-01',
                                       tz = 'UTC')
      
      #match the layer time to the nearest hour.
      #' `inner.layer` first.
      inner.layer.times <- gsub('X', '', names(inner.emissions.raster))
      inner.layer.times <- as.POSIXlt(as.numeric(inner.layer.times),
                                      origin = '1970-01-01',
                                      tz = 'UTC')
      min.idx <- which.min(inner.layer.times - tmp.foot.layer.pos)
      
      #grab the inner domain's hourly emissions
      tmp.inner.layer <- inner.emissions.raster[[min.idx]]
      
      #' `outer.layer` next.
      outer.layer.times <- gsub('X', '', names(outer.emissions.raster))
      outer.layer.times <- as.POSIXlt(as.numeric(outer.layer.times),
                                      origin = '1970-01-01',
                                      tz = 'UTC')
      min.idx <- which.min(outer.layer.times - tmp.foot.layer.pos)
      
      #grab the outer domain's hourly emissions
      tmp.outer.layer <- outer.emissions.raster[[min.idx]]
      
      #same as above with bio layers
      if(!is.na(biosphere.emissions)) {
        
        bio.layer.times <- gsub('X', '', names(bio.emissions.raster))
        bio.layer.times <- as.POSIXlt(as.numeric(bio.layer.times),
                                      origin = '1970-01-01',
                                      tz = 'UTC')
        min.idx <- which.min(bio.layer.times - tmp.foot.layer.pos)
        
        #grab the biospheric hourly emissions
        tmp.bio.layer <- bio.emissions.raster[[min.idx]]

      }
      
      #determine the domain's total emissions
      tmp.total.layer <-
        extend(tmp.inner.layer, tmp.outer.layer, value = 0) +
        tmp.outer.layer
      
      #add the (optional) biospheric layer
      if(!is.na(biosphere.emissions))
        tmp.total.layer <- tmp.total.layer + tmp.bio.layer
      
      #inner domain emissions
      hourly.XCO2[j] <-
        suppressWarnings(cellStats(tmp.footprint.timestep*tmp.total.layer,
                                   sum))
      
      #outer domain emissions
      hourly.bkg.XCO2[j] <-
        suppressWarnings(cellStats(tmp.footprint.timestep*tmp.outer.layer,
                                   sum))
      
    }; XCO2 <- sum(hourly.XCO2); bkg.XCO2 <- sum(hourly.bkg.XCO2)
    
    #include the error here
    #how the error is added depends on how it's quantified
    if(!is.na(obs.error)) {
      
      #' If an observation is provided in the main script,
      #' *inversion.analysis.r* then the other error method
      #' will be overridden.
      collect.error <- rnorm(1, 0, obs.error)
      
    } else if(is.na(obs.error)) {
      
      #mean and s.d. for a normal distribution
      if(!all(is.na(errors[,2])) & !all(is.na(errors[,3]))) {
        
        #add standard deviations in quadrature
        combined.sd <- sqrt(sum(errors[,3]^2))
        
        # Obtain the error value from a normal distribution
        collect.error <- rnorm(1, 0, combined.sd)
        
        #max percentage of value used in a uniform distribution
      } else if((all(is.na(errors[,2])) & all(is.na(errors[,3]))) &
                all(is.na(errors[,5])) & !all(is.na(errors[j,4]))) {
        
        #add percent errors in quadrature
        combined.percent <- sqrt(sum(errors[,4]^2))
        
        #Obtain the error value from a percent error
        collect.error <- runif(1,
                               -1*combined.percent*XCO2,
                               combined.percent*XCO2)
        
        #RMSE used as s.d. in normal distribution (mu = 0)
      } else if((all(is.na(errors[,2])) & all(is.na(errors[,3])) &
                 all(is.na(errors[,4]))) & !all(is.na(errors[,5]))) {
        
        #add the RMSE values in quadrature
        combined.rmse <- sqrt(sum(errors[,5]^2))
        
        #start adding the error term for a normal distribution
        collect.error <- rnorm(1, 0, combined.rmse)
        
      } else{stop('Check error file!')}
    }
    
    ### Inner domain - obs.rds
    #add to the XCO2 dataframe with error
    add.line <-
      data.frame(tmp.date.time,
                 as.numeric(tmp.lon),
                 as.numeric(tmp.lat),
                 XCO2 + collect.error)
    names(add.line) <- names(XCO2_df)
    
    #add the line to the dataframe
    XCO2_df <- rbind(add.line, XCO2_df)
    
    ### Background domain - background.rds
    add.line <-
      data.frame(tmp.date.time,
                 as.numeric(tmp.lon),
                 as.numeric(tmp.lat),
                 bkg.XCO2)
    names(add.line) <- names(XCO2_df)
    
    #add the line to the dataframe (background)
    bkg.XCO2_df <- rbind(add.line, bkg.XCO2_df)
    
  } #closes the footprint loop
  
  #generate obs.rds
  obs.rds <-
    as.matrix(data.frame(site_code = site_code.list,
                         seconds_since_1970_01_01 =
                           paste(as.numeric(start.time)),
                         co2_ppm = paste(XCO2_df$xco2)))
  saveRDS(obs.rds, file = obs.output)
  
  if(std.bkg == F) {
    
    stop('Non-standard background not currently supported!')
    
    # ### generate background.rds ###
    # #convert to epoch time
    # epoch.time <-
    #   as.POSIXlt(bkg.XCO2_df$date.time,
    #              format = '%Y%m%d%H%M',
    #              origin = '1970-01-01',
    #              tz = 'UTC')
    # bkg.XCO2_df$date.time <- epoch.time
    # names(bkg.XCO2_df) <- c('time', 'lon', 'lat', 'avg_CO2')
    # saveRDS(bkg.XCO2_df, file = background.output)
    
  } else if(std.bkg == T) {
    
    #sets background to zero
    bkg.XCO2_df <- data.frame(time = start.time,
                              avg_CO2 = 0)
    saveRDS(bkg.XCO2_df, file = background.output)
  }
  
}
