# subroutine to readin tiff format of 1kmx1km ODIAC emissions
# and return *.nc format with selected emissions given a lat/lon domain
# This is a shameless knockoff of D. Wu's function tif2nc.odiacv3.r
# DW, update with ODIACv2017, 11/01/2017

get.odiac <- function(tiff.path, nc.extent, YYYYMM = NULL,
                      convert.units = TRUE) {

  # Some QA/QC
  if(!is.character(YYYYMM) & !is.null(YYYYMM))
    stop('YYYYMM must be a STRING!')
  
  # Load required libraries
  library(raster)
  
  YYYY <- as.numeric(substr(YYYYMM, 1, 4))
  YY_s <- substr(YYYYMM, 3, 4) #string format for searching
  MM <- as.numeric(substr(YYYYMM, 5, 6))
  MM_s <- substr(YYYYMM, 5, 6) #string format for searching
  
  # Find the closest match
  odiac.file <-
    list.files(tiff.path, pattern = paste0(YY_s, MM_s), full.names = TRUE)
  if(length(odiac.file) == 0) {
    message('No exact match for ODIAC file. Matching by month instead.')
    odiac.file <- list.files(tiff.path,
                             pattern = paste0(MM_s, '.tif'),
                             full.names = TRUE)
  }
  
  if(length(odiac.file) == 0)
    message('Cannot find appropriate ODIAC file.')
  
  # 21600 rows and 43200 columns
  emiss <- raster(odiac.file) # convert to raster
  
  # Check to see if YYYY is a leap year
  # Difference in second entry of each vector (February)
  if((YYYY %% 400 == 0) | (YYYY %% 4 == 0)) {
    month.days <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[MM]
  } else {
    month.days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[MM]
  }

  # subset spatial domain
  sel.emiss <- crop(emiss, nc.extent)
  remove('emiss') #save RAM

  if(convert.units) {
  
    # Method 2 -- compute area using area() function in raster package
    area.raster <- raster::area(sel.emiss) * 1E6    # convert km2 to m2
    
    # convert the unit of CO2 emiss from Tonne Carbon/cell/month to umol/m2/s
    sel.emiss <- sel.emiss * 1E6 / 12 * 1E6 # convert tonne-C to uomol-C (= umole-CO2)
    sel.emiss <- sel.emiss / month.days / 24 / 60 / 60	# convert per month to per second
    sel.emiss <- sel.emiss / area.raster		# convert per cell to per m2
    # NOW sel.co2 has unit of umole-CO2/m2/s, can be used directly with footprint
    
  }

  # finally, return raster
  return(sel.emiss)
}
