timestep.normalization.constant <- function(lon.lat, sector_df,
                                            ODIAC,
                                            carma.inventory,
                                            edgar.inventory, nc.extent,
                                            temporal.downscaling,
                                            time) {

  #Be sure that ODIAC does not have converted units!  

  #break the ODIAC file into non-LPS file
  ODIAC.nlps <- get.odiac.nlps(CarMA = read.csv(carma.inventory), ODIAC)
  remove('ODIAC') #save RAM
  
  #category loop
  for(i in 1:nrow(sector_df)) {
    
    cat.name <- sector_df[i,1]
    cat.sectors <- unlist(str_split(sector_df[i,2], pattern = ' '))
    
    for(j in 1:length(cat.sectors)) {
      
      #large point sources (LPSs) first.
      #DON'T use them in the normalization constant!
      #These are dealt with differently.
      
      #construct the sector path and grab the static EDGAR sector
      sector.path <- list.files(edgar.inventory,
                                pattern = cat.sectors[j],
                                full.names = TRUE)
      
      edgar.sector <- get.edgar.sector(nc.path = sector.path,
                                       nc.extent = extent(ODIAC.nlps),
                                       annual.total = FALSE)
      resampled.edgar.sector <- resample(edgar.sector, ODIAC.nlps)
      resampled.edgar.sector[is.na(resampled.edgar.sector)] <- 0
      resampled.edgar.sector[resampled.edgar.sector < 0] <- 0
      remove('edgar.sector') #save RAM
      
      #get the sector weighting
      weight <- edgar.sector.weighting(citylon = lon.lat$citylon,
                                       citylat = lon.lat$citylat,
                                       local.tz = lon.lat$tz,
                                       sector.name = cat.sectors[j],
                                       temporal.downscaling.files =
                                         temporal.downscaling,
                                       time = time,
                                       monthly = TRUE)
      
      #construct with the LPS version of ODIAC
      weighted.edgar.sector <- weight*resampled.edgar.sector
      remove('resampled.edgar.sector') #save RAM
      
      # This EDGAR raster is in kg / m2 / hr.
      # convert the unit of CO2 emiss from kg/m2/hr to umol/m2/s
      weighted.edgar.sector <- weighted.edgar.sector * (10^9 / 44.01) # convert kg to uomol-C (= umole-CO2)
      weighted.edgar.sector <- weighted.edgar.sector / 60 / 60	# convert per hour to per second
      # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
      
      #Large point sources are handled differently.
      if(cat.sectors[j] == 'ENE')
        weighted.edgar.sector <- 0*weighted.edgar.sector
      
      if(j == 1)
        category.flux.layer <- weighted.edgar.sector
      
      if(j > 1)
        category.flux.layer <-
        category.flux.layer + weighted.edgar.sector
      
    } #closes the sectors making up each category
    
    if(i == 1) total.flux <- category.flux.layer
    if(i > 1) total.flux <- total.flux + category.flux.layer
    
  } #closes categories
  
  return(total.flux)
  
} #closes function