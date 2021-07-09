edgar.sector.weighting <- function(citylon = NULL, citylat = NULL, local.tz = NULL,
                                   sector.name = NULL, temporal.downscaling.files = NULL,
                                   time = NULL, monthly = FALSE) {
  
  # First, read in the butt load of temporal downscaling csv files
  scaling.files <- list.files(temporal.downscaling.files, full.names = TRUE,
                              pattern = '\\.csv$')
  
  if(length(scaling.files) != 8)
    stop('Incorrect number of temporal downscaling files.')
  
  # Read in each csv file and name give it a variable name based on
  # the csv base name. There should be 7 csv files for EDGARv5.
  for(i in 1:8) {
    scaling.file.name <- substr(basename(scaling.files[i]), 1,
                                nchar(basename(scaling.files[i]))-4)
    scaling.file.path <- scaling.files[i]
    eval(parse(text = paste0(scaling.file.name,
                             " <- read.csv('", scaling.file.path, "', sep = ';')")))
  }
  
  # Before running through the temporal downscaling files,
  # grab the parent activity code.
  activity_code <- unlist(strsplit(sector.name, split = '_'))[1]
  
  #' First, determine the time zone (tz) code
  #' Data from `lat_lon_TZ_id.csv`
  tz.idx <- which.min(sqrt((lat_lon_TZ_id$lon - citylon)^2 +
                             (lat_lon_TZ_id$lat - citylat)^2))
  TZ_ID <- lat_lon_TZ_id$TZ_ID[tz.idx]
  
  #' Next, identify the country code and UTC_reference.
  #' Data from `timezones_definition.csv`
  Country_code_A3 <- subset(timezones_definition, TZ_id == TZ_ID)$Country_code_A3
  UTC_reference <- subset(timezones_definition, TZ_id == TZ_ID)$UTC_reference
  
  #' The weekend type is needed. This changes depending on the country
  #' Data from `weekenddays.csv`
  cntry <- Country_code_A3 #rename for subset
  Weekend_type_id <- subset(weekenddays, Country_code_A3 == cntry)$Weekend_type_id
  
  #' Identify the Daytype_id of each day of the assigned week
  #' Data from `weekdays.csv`
  wk_typ_id <- Weekend_type_id #rename for subset
  sub.weekdays <- subset(weekdays, Weekend_type_id == wk_typ_id)
  
  #' Identify the daily_factor of each day of the assigned week
  #' Data from `weekly_profiles.csv`
  
  # as.POSIXlt iterates Sun thru Sat as 0 thrus 6. Sun must become 7 for EDGAR
  # Convert the UTC time to local time
  POSIX.format <- as.POSIXlt(time, format = '%Y.%m.%d.%H.%M.%S', tz = 'UTC')
  UTC_shift <- suppressWarnings(tz_offset(POSIX.format, tz = local.tz))$utc_offset_h
  native.POSIX.format <- POSIX.format + UTC_shift*3600
  
  # Get the weekday name and number
  Weekday_name <- weekdays(as.Date(native.POSIX.format))
  Weekday_id <- subset(sub.weekdays, weekday_name == Weekday_name)$Weekday_id
  
  # Monthly profiles are not currently provided for EDGAR v5.
  # A simple division by 12 is used for now.
  MONTHLY_FACTOR <- 1/12

  # For some unknown reason, DAILY_FACTOR for sector PRO has multiple entries.
  act.cde <- activity_code; Wkdy_id <- Weekday_id #rename for subset
  DAILY_FACTOR <- subset(weekly_profiles,
                         Country_code_A3 == cntry &
                           activity_code == act.cde &
                           Weekday_id == Wkdy_id)$daily_factor[1]
  if(length(DAILY_FACTOR) == 0 | is.na(DAILY_FACTOR)) DAILY_FACTOR <- 0
  
  #' Determine the Daytype_id by using the Weekend_id (`Wkdy_id`)
  Dy_typ <- subset(sub.weekdays, Weekday_id == Wkdy_id)$Daytype_id
  
  #' Identify the hourly_factor of each day of the assigned week
  #' Data from `hourly_profiles.csv`
  #' For now, the factor from the nearest hour is applied.
  #' This needs to be linearly interpolated later, with particular
  #' attention given to the end times of each day.
  hr <- as.POSIXlt(native.POSIX.format)$hour
  min <- as.POSIXlt(native.POSIX.format)$min
  
  # Linear interpolation of the hourly factor
  # initial time
  HOURLY_FACTOR_1 <- subset(hourly_profiles,
                          Country_code_A3 == cntry &
                            activity_code == act.cde &
                            month_id == month(POSIX.format) &
                            Daytype_id == Dy_typ)[,(4+hr)]
  
  # final time
  HOURLY_FACTOR_2 <- subset(hourly_profiles,
                            Country_code_A3 == cntry &
                              activity_code == act.cde &
                              month_id == month(POSIX.format) &
                              Daytype_id == Dy_typ)[,(5+hr)]
  
  if(length(HOURLY_FACTOR_1) == 0 | length(HOURLY_FACTOR_2) == 0) {
    HOURLY_FACTOR <- 0
  } else {
    # interpolate here
    HOURLY_FACTOR <-
      HOURLY_FACTOR_1 + (min/60)*(HOURLY_FACTOR_2 - HOURLY_FACTOR_1)
  }
  
  #' Apply the temporal downscaling to the EDGAR sector here
  #' Determine days in month for part of the downscaling process
  #' Resample the sector based on the supplied footprint (bilinear interp.)
  #' After the bilinear interpolation, replace negative values with zero.
  #' Weighted as described in `Crippa et al. (2020)`.
  n_days <- days_in_month(POSIX.format)
  
  if(!monthly) {
    weighting <- (MONTHLY_FACTOR*(7/n_days)*DAILY_FACTOR*HOURLY_FACTOR)
    names(weighting) <- 'weight.value'
  } else if(monthly) {
    weighting <- (7/n_days)*DAILY_FACTOR*HOURLY_FACTOR
    names(weighting) <- 'weight.value'
  }
  
  return(weighting)
  
}