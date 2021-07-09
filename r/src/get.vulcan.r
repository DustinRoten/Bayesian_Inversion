get.vulcan <- function(vulcan.path = NULL, sector = 'total',
                       year = NA, loc = 'US', ref.raster = NA) {
  
  loc.list <- list.files(vulcan.path, pattern = sector, full.names = TRUE)
  
  loc.idx <- which(!is.na(str_match(loc.list, pattern = loc)))
  
  vulcan.file <-
    list.files(vulcan.path, pattern = sector, full.names = TRUE)[loc.idx]
  vulcan.raster <- stack(vulcan.file)
  
  #get the appropriate layer of Vulcan
  list.years <- names(vulcan.raster)
  years.as.numeric <- as.numeric(substr(list.years, 2, 5))
  
  #find the closest year
  which.layer <- which.min(abs(years.as.numeric - year))
  layer.name <- list.years[which.layer]
  eval(parse(text = paste0('Vulcan <- vulcan.raster$', layer.name)))
  remove('vulcan.raster')
  
  if(length(ref.raster) == 1) {
    
    message(paste0('Returning original Vulcan layer ',
                   list.years[which.layer]))
    return(Vulcan)
    
  } else {
    
    new.Vulcan <- projectRaster(Vulcan, ODIAC)
    remove('Vulcan')
    
    # Check to see if Vulcan year is a leap year
    # Difference in second entry of each vector (February)
    if((years.as.numeric[which.layer] %% 400 == 0) |
       (years.as.numeric[which.layer] %% 4 == 0)) {
      year.seconds <- 366*24*60*60
    } else {
      year.seconds <- 365*24*60*60
    }
    
    #convert units (code from Dien Wu)
    # convert the unit of CO2 emiss from tC/km2/year to umol/m2/s
    new.Vulcan <- new.Vulcan*1E6/12*1E6 # convert tonne-C to umol-C (= umole-CO2)
    new.Vulcan <- new.Vulcan/year.seconds # convert per year to per second
    new.Vulcan <- new.Vulcan/1E6	# convert per cell to per m2
    # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
    
    return(new.Vulcan)
  }
  
}
