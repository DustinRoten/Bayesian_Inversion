get_vulcan_emiss_hourly <- function(vulcan_emiss = NULL, category = NA,
                                    date.time = NA, loc = 'US', ref.raster = NA,
                                    include.maritime = TRUE, include.elec_prod = TRUE) {
  
  require(stringr)
  
  if(is.na(category)) {
    
    if(include.maritime & include.elec_prod) {
      #list all available rasters
      all.rasters <- list.files(vulcan_emiss, pattern = 'total',
                                full.names = TRUE)
    } else if(!include.maritime | !include.elec_prod) {
      #' If `include.maritime == F`, then it makes the most
      #' sense to subtract the maritime sector (*cmv*);
      #' however, in an effort to keep the code similar to 
      #' other sections, all sectors except the maritime
      #' sector will be summed.
      if(!include.maritime & include.elec_prod)
      all.rasters <- grep(list.files(vulcan_emiss, full.names = TRUE),
                          pattern = 'total|cmv', invert = TRUE,
                          value = TRUE)
      
      #' If `include.elec_prod == F`, then it makes the most
      #' sense to subtract the electricity production sector;
      #' (*elec_prod*) however, in an effort to keep the code
      #' similar to other sections, all sectors except the
      #' electricity production sector will be summed.
      if(include.maritime & !include.elec_prod)
        all.rasters <- grep(list.files(vulcan_emiss, full.names = TRUE),
                            pattern = 'total|elec_prod', invert = TRUE,
                            value = TRUE)
      
      #' A combination of the two cases above
      if(!include.maritime & !include.elec_prod)
        all.rasters <- grep(list.files(vulcan_emiss, full.names = TRUE),
                            pattern = 'total|elec_prod|cmv',
                            invert = TRUE, value = TRUE)
    }
    
    #subset by locations
    all.rasters.idx <-
      which(!is.na(str_match(all.rasters, pattern = loc)))
    all.rasters <- all.rasters[all.rasters.idx]
    
    #required year
    year.number <- year(date.time)
    
    #From the identified categories, get the available file paths
    idx <- grep(pattern = paste0(year.number),
                x = basename(all.rasters))
    available.rasters <- all.rasters[idx]
    if(length(available.rasters) == 0) {
      i <- 0
      while(length(available.rasters) == 0) {
        i <- i+1
        idx <- grep(pattern = paste0(year.number - i),
                    x = basename(all.rasters))
        available.rasters <- all.rasters[idx]
      }
      cat(paste0('Using Vulcan year ', year.number - i, ' instead.'),
          '\n')
    }; remove('i')
    
    #determine the day number in the year
    #(there is probably a fancier way to do this...)
    day.number <- as.POSIXlt(date.time)$yday
    if(nchar(day.number) == 1) day.string <- paste0('d00', day.number)
    if(nchar(day.number) == 2) day.string <- paste0('d0', day.number)
    if(nchar(day.number) == 3) day.string <- paste0('d', day.number)
    
    #match the file path with the required day
    which.raster <-
      str_match(available.rasters,
                pattern = day.string)
    which.raster <- which(!is.na(which.raster))
    
    #retrieve each raster and layer names
    for(i in 1:length(which.raster)) {
      
      vulcan.raster <- brick(available.rasters[which.raster[i]])
      
      pos.layer.names <- as.POSIXlt(gsub('X', '', names(vulcan.raster)),
                                    format = '%Y.%m.%d.%H.%M.%S',
                                    tz = 'UTC')
      
      #match the hours to detemrine which layer to use
      hour.idx <- which(hour(pos.layer.names) == hour(date.time))
      
      #grab the layer and add it to the stack
      if(i == 1) {
        tmp.layer <- vulcan.raster[[hour.idx]]
        tmp.layer[is.na(tmp.layer)] <- 0 #remove any NA's in the sectors
        vulcan.layer <- tmp.layer
      } else if(i > 1) {
        tmp.layer <- vulcan.raster[[hour.idx]]
        tmp.layer[is.na(tmp.layer)] <- 0 #remove any NA's in the sectors
        vulcan.layer <- vulcan.layer + tmp.layer
      }
    
    }; remove('vulcan.raster') #conserve RAM
    
  } else if(!is.na(category)) {
    
    stop('Sector-based Vulcan is not yet available!')
    
    #############################################################################
    #' Here, a matching layer must be added for each timestep in the prior_emiss.
    #' Don't forget to add conditional statements to deal with each emissions
    #' category!!! (A separate CSV file will be required.)
    #############################################################################
    
  } #ends categories
  
  #Convert units and return the vulcan layer
  if(class(ref.raster) == 'logical') {
    
    message(paste0('Returning original Vulcan layer.'))
    return(vulcan.layer)
    
  } else {
    
    #Reproject the Vulcan raster onto the ref.raster (if applicable)
    if(class(ref.raster) == "RasterLayer" | class(ref.raster) == "RasterBrick")
      Vulcan <- projectRaster(vulcan.layer, ref.raster)
    
    ### convert units (code from Dien Wu)
    #convert the unit of CO2 emiss from tonne-C/km2/hr to umol/m2/s
    Vulcan <- Vulcan*(1E6/12*1E6) #First, convert tonne-C to umol-C (= umole-CO2)
    Vulcan <- Vulcan/3600 # convert per hour to per second
    Vulcan <- Vulcan/1E6	# convert per cell to per m2
    # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
    
  }
  
  return(Vulcan)
}