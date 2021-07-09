get.edgar.total <- function(edgar.inventory, downscaling.extension,
                            lon.lat, categories, time, monthly = TRUE,
                            ODIAC, convert.units = TRUE) {
  
  #list all of the sectors required for this study
  all.sectors <- NULL
  for(i in 1:nrow(categories)) {
    cat.sectors <- unlist(strsplit(categories[i,2], split = ' '))
    all.sectors <- c(all.sectors, cat.sectors)
  }
  
  for(i in 1:length(all.sectors)) {
    
    #identify the sector name
    sector.name <- all.sectors[i]
    
    #construct the sector path
    edgar.sector.path <- list.files(edgar.inventory,
                                    pattern = sector.name,
                                    full.names = TRUE)
    
    #get the edgar sector
    edgar.sector <- get.edgar.sector(nc.path = edgar.sector.path,
                                     nc.extent = extent(ODIAC))
    edgar.sector <- resample(edgar.sector, ODIAC)
    edgar.sector[edgar.sector < 0] <- 0
    
    #calculate the weighting
    weight <- edgar.sector.weighting(citylon = lon.lat$citylon,
                                     citylat = lon.lat$citylat,
                                     local.tz = lon.lat$tz,
                                     sector.name = all.sectors[i],
                                     temporal.downscaling.files =
                                       downscaling.extension,
                                     time, monthly)
    
    #calculate the weighted raster
    weighted.edgar.sector <- weight*edgar.sector
    remove('edgar.sector') #save RAM
    
    if(i == 1) total.edgar <- weighted.edgar.sector
    if(i > 1) total.edgar <- total.edgar + weighted.edgar.sector
    
  } #close all.sectors loop
    
  if(convert.units) {
    
    ### Convert the weighted layer to umol/m2/s ###
    # Compute area using area() function in raster package
    area.raster <- raster::area(total.edgar) * 1E6 # convert km2 to m2
    
    # Convert the unit of CO2 emiss from Tonne Carbon/cell/hour to umol/m2/s
    total.edgar <- total.edgar * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
    total.edgar <- total.edgar / 60 / 60	# convert to per second
    total.edgar <- total.edgar / area.raster		# convert per cell to per m2
    # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
    
    remove('area.raster') #save RAM
    
  }
  
  return(total.edgar)
  
}
  