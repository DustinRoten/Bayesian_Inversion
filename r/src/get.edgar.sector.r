get.edgar.sector <- function(nc.path = NULL, nc.extent = NULL,
                             annual.total = FALSE, year = NA) {
  
  library(lubridate)
  
  if(is.na(year)) year <- as.numeric(basename(dirname(nc.path)))
  
  if(!annual.total) {
    
    # Read in the original EDGAR sector (annual)
    # Convert to umol/m2/s instead of kg/m2/s
    original.edgar.sector <- (10^9/44.01)*raster(nc.path)
    
    # Break the raster into halves and rearrange
    # Merge the two halves back together and update extents
    # Crop the result based on the supplied extent (nc.extent)
    left.edgar <- crop(original.edgar.sector, c(0, 180, -90, 90))
    right.edgar <- crop(original.edgar.sector, c(180, 360, -90, 90))
    extent(right.edgar) <- c(-180, 0, -90, 90)
    edgar.sector <- merge(right.edgar, left.edgar)
    edgar.sector <- crop(edgar.sector, nc.extent)
    
  } else if(annual.total) {
    
    # Leave in kg/m2/s
    original.edgar.sector <- raster(nc.path)
    
    # Break the raster into halves and rearrange
    # Merge the two halves back together and update extents
    # Crop the result based on the supplied extent (nc.extent)
    left.edgar <- crop(original.edgar.sector, c(0, 180, -90, 90))
    right.edgar <- crop(original.edgar.sector, c(180, 360, -90, 90))
    extent(right.edgar) <- c(-180, 0, -90, 90)
    edgar.sector <- merge(right.edgar, left.edgar)
    edgar.sector <- crop(edgar.sector, nc.extent)
    
    #convert to annual amount
    is.leap.year <- leap_year(year)
    if(is.leap.year) {days <- 366} else {days <- 365}
    
    #seconds in a year
    seconds <- days*24*3600
    edgar.sector <- edgar.sector*seconds
  }
  
  return(edgar.sector)
  
}