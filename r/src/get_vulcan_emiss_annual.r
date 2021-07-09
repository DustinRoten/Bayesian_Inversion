get_vulcan_emiss_annual <- function(vulcan_emiss = NULL, category = NA,
                                    years = NA, loc = 'US', ref.raster = NA) {
  
  require(stringr)
  
  if(is.na(category)) {
  
    #get a list of total (annual) Vulcan files
    loc.list <-
      list.files(vulcan_emiss, pattern = 'total', full.names = TRUE)
    contains.us <- str_match(loc.list, pattern = 'US')
    loc.idx <- which(!is.na(contains.us))
    
    #read in the vulcan raster
    vulcan.raster <- brick(loc.list[loc.idx])
    
    #####
    #get the appropriate layer of Vulcan by matching the years
    list.years <- names(vulcan.raster)
    years.as.numeric <- as.numeric(substr(list.years, 2, 5))
    
    #find the closest year for each prior_emiss year
    which.layers <- NULL
    for(i in 1:length(years)) {
      which.layers[i] <- which.min(abs(years.as.numeric - years[i]))
    }
    #' Some layers may have been duplicated in the above loop.
    #' Only require unique layers, no duplicates.
    layer.name <- list.years[unique(which.layers)]
    eval(parse(text = paste0('Vulcan <- vulcan.raster$', layer.name)))
    remove('vulcan.raster')
    
    #Reproject the Vulcan raster onto the ref.raster (if applicable)
    if(class(ref.raster) == "RasterLayer" | class(ref.raster) == "RasterBrick")
      Vulcan <- projectRaster(Vulcan, ref.raster)
    #####
  
  } else if(!is.na(category)) {
    
    stop('Sector-based Vulcan is not yet available!')
    
    #############################################################################
    #' Here, a matching layer must be added for each timestep in the prior_emiss.
    #' Don't forget to add conditional statements to deal with each emissions
    #' category!!! (A separate CSV file will be required.)
    #############################################################################
    
  }
  
  if(class(ref.raster) == 'logical') {
    
    message(paste0('Returning original Vulcan layer ',
                   list.years[which.layer]))
    return(Vulcan)
    
  } else {
    
    #setup an output raster to attach converted layers to.
    #add a dummy layer of all zeros. This will be removed later.
    output.raster <- Vulcan[[1]]
    values(output.raster) <- 0
    for(i in 1:length(which.layers)) {
      
      #step through each layer of Vulcan
      tmp.layer <- which.layers[i]
      
      # Check to see if Vulcan year is a leap year
      # Difference in second entry of each vector (February)
      if((years.as.numeric[tmp.layer] %% 400 == 0) |
         (years.as.numeric[tmp.layer] %% 4 == 0)) {
        year.seconds <- 366*24*60*60
      } else {
        year.seconds <- 365*24*60*60
      }
    
      ### convert units (code from Dien Wu)
      #convert the unit of CO2 emiss from tC/km2/year to umol/m2/s
      #first, get the layer name
      convert.layer <- paste0('Vulcan$', list.years[tmp.layer])
      
      # convert tonne-C to umol-C (= umole-CO2)
      eval(parse(text = paste0('tmp.raster.layer <- ',
                               convert.layer, '*1E6/12*1E6')))
      
      tmp.raster.layer <- tmp.raster.layer/year.seconds # convert per year to per second
      tmp.raster.layer <- tmp.raster.layer/1E6	# convert per cell to per m2
      # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
      
      output.raster <- addLayer(output.raster, tmp.raster.layer)
     
    }
    
    #remove the dummy layer and update the names here
    output.raster <- dropLayer(output.raster, i = 1)
    names(output.raster) <- list.years[which.layers]
    
    return(output.raster)
  }
  
}
