setup.config <- function(timesteps = NULL, src_path = NULL, out_path = NULL, foot_dir = NULL,
                         include_dir = NULL, clear_H = TRUE, compute_chi_sq = FALSE,
                         include_outer = TRUE, include_bio = TRUE, lonlat_domain_file.NAME = NULL,
                         lonlat_outer_file.NAME = NULL, prior_file.NAME = NULL, sigma_file.NAME = NULL,
                         outer_file.NAME = NULL, bio_file.NAME = NULL, obs_file.NAME = NULL,
                         bg_file.NAME = NULL, lon_res = NULL, lat_res = NULL, foot_dim_lonxlat = T,
                         aggregate_obs = FALSE, min_agg_obs = 2, ls = 1, lt = 1, flux_units = "umol/(m2 s)",
                         save.to = NULL) {
  
  #change the working directory
  setwd(save.to)
  
  # ~~~~~~~~~~~~~~ Necessary paths ~~~~~~~~~~~~~~#
  
  # Names of observational sites
  #(correspond to directory names in the footprint dir)
  sites <- list.files(foot_dir) #modified by DR 6/1/2021
  nsites <- length(sites)
  
  # ~~~~~~~~~~~~~~ Info files ~~~~~~~~~~~~~~#
  
  # lonlat_domain file - lists all lon/lat pairs in domain as look up table
  lonlat_domain_file <- paste0(include_dir, lonlat_domain_file.NAME)
  
  # lonlat_outer file - lists all lon/lat pairs in outer domain as look up table
  # NOTE: can be NA if include_outer and include_bio are both FALSE
  lonlat_outer_file <- paste0(include_dir, lonlat_outer_file.NAME)
  
  # Prior emissions file
  prior_file <- paste0(include_dir, prior_file.NAME)
  
  # Prior uncertainty file
  sigma_file <- paste0(include_dir, sigma_file.NAME)
  
  # outer-domain emissions inventory file
  # NOTE: can be NA if include_outer is FALSE
  outer_file <- paste0(include_dir, outer_file.NAME)
  
  # Biospheric flux inventory file
  # NOTE: can be NA if include_bio is FALSE
  bio_file <- paste0(include_dir, bio_file.NAME)
  
  
  # obs file should contain data for all sites used in the inversion, with
  # nomenclature:
  # SITECODE  OBS TIME (seconds since 1970-01-01 00:00:00)  OBS VALUE
  # 'site1'          '[t1 number of seconds]'                 'X.XXX'
  # 'site1'          '[t2 number of seconds]'                 'X.XXX'
  #   ...                       ...                             ...
  # 'site2'          '[t1 number of seconds]'                 'X.XXX'
  
  obs_file <- paste0(include_dir, obs_file.NAME)
  
  # bg filepath (follows same convention as obs, but without column 1)
  bg_file <- paste0(include_dir, bg_file.NAME)
  
  
  # ~~~~~~~~~~~~~~ Model params ~~~~~~~~~~~~~~#
  
  # spatial grid resolution
  round_digs <- 1 + max(nchar(strsplit(as.character(lon_res), ".", fixed = TRUE)[[1]][[2]]),
                        nchar(strsplit(as.character(lat_res), ".", fixed = TRUE)[[1]][[2]]))
  
  # flux time range
  flux_year_start <- year(min(timesteps))
  flux_month_start <- month(min(timesteps))
  flux_day_start <- day(min(timesteps))
  flux_hour_start <- hour(min(timesteps))
  flux_min_start <- minute(min(timesteps))
  
  # flux bins reflect [start -> end) i.e. flux end time marks the end of last bin,
  # not start of last bin (unlike obs)
  flux_year_end <- year(max(timesteps))
  flux_month_end <- month(max(timesteps))
  flux_day_end <- day(max(timesteps))
  flux_hour_end <- hour(max(timesteps))
  flux_min_end <- minute(max(timesteps))
  
  flux_t_res <- abs(min(diff(as.numeric(timesteps))))  #flux time resolution, in seconds
  
  ### Added by D. Roten 6/3/2021 ###
  #read in the observation file to set up the following parameters
  obs.data <- readRDS(file.path(save.to, obs_file))
  obs.times <- as.POSIXlt(as.numeric(obs.data[,2]),
                          origin = '1970-01-01',
                          tz = 'UTC')
  obs.times <- unique(obs.times)
  
  # obs/receptor time range
  obs_year_start <- year(min(obs.times))
  obs_month_start <- month(min(obs.times))
  obs_day_start <- day(min(obs.times))
  obs_hour_start <- hour(min(obs.times))
  obs_min_start <- minute(min(obs.times))
  
  # obs bins reflect [start -> end] i.e. obs end time marks beginning of last bin.
  obs_year_end <- year(max(obs.times))
  obs_month_end <- month(max(obs.times))
  obs_day_end <- day(max(obs.times))
  obs_hour_end <- hour(max(obs.times))
  obs_min_end <- minute(max(obs.times))
  
  #observation time resolution, in seconds
  obs_t_res <-
    suppressWarnings(abs(min(diff(as.numeric(obs.times)))))
  if(length(obs.times) == 1) obs_t_res <- 0
  
  # if observations are only desired from a subset of time, declare here otherwise,
  # set to 0:23 (i.e. all hours of day)
  subset_hours_utc <- 0:23
  
  # get the number of days
  ndays <- length(seq(from = ISOdatetime(obs_year_start, obs_month_start, obs_day_start,
                                         obs_hour_start, obs_min_start, 0, tz = "UTC"),
                      to = ISOdatetime(obs_year_end, obs_month_end, obs_day_end,
                                       obs_hour_end, obs_min_end, 0, tz = "UTC"),
                      by = 24 * 3600))
  
  # Compute the number of time steps
  ntimes <- length(seq(from = ISOdatetime(flux_year_start, flux_month_start, flux_day_start,
                                          flux_hour_start, flux_min_start, 0, tz = "UTC"),
                       to = ISOdatetime(flux_year_end, flux_month_end, flux_day_end,
                                        flux_hour_end, flux_min_end, 0, tz = "UTC"),
                       by = flux_t_res)) # -1 removed by D. Roten 6/3/21
  
  # time steps for posterior covariance V_shat, the timesteps to aggregate over
  tstart <- 1
  tstop <- ntimes
  
  # define string for flux units
  flux_units <- "umol/(m2 s)"
  
  # ~~~~~~~~~~~ Clear all existing data files before proceeding ~~~~~~~~~#
  
  # get all files in the outer directory
  out_files <- paste0(out_path, list.files(out_path))
  
  if(!dir.exists("H/"))
    dir.create("H/")
  
  if(!is.na(lonlat_outer_file) & !dir.exists("H_outer/"))
    dir.create("H_outer/")
  
  if(!dir.exists("HQ/"))
    dir.create("HQ/")
  
  if(!dir.exists("out/"))
    dir.create("out/")
  
  
  if (clear_H) {
    # H files
    if (length(list.files("H/")) > 0) {
      out_files <- c(out_files, paste0("H/", list.files("H/")))
    }
    
    # H files (outer domain)
    if (!is.na(lonlat_outer_file) & length(list.files("H_outer/")) > 0) {
      out_files <- c(out_files, paste0("H_outer/", list.files("H/")))
    }
  }
  
  # HQ files
  if (length(list.files("HQ/")) > 0) {
    out_files <- c(out_files, paste0("HQ/", list.files("HQ/")))
  }
  
  # determine which of these files exist (redundant but good)
  iFiles <- which(file.exists(out_files))
  
  # remove existing files to clean the directory
  invisible(sapply(out_files[iFiles], FUN = function(x) system(paste("rm", x))))
  
  save(src_path, out_path, foot_dir, include_dir,
       sites, nsites,
       clear_H, compute_chi_sq, include_outer, include_bio,
       lonlat_domain_file, lonlat_outer_file, prior_file,
       sigma_file, outer_file, bio_file, obs_file, bg_file,
       lon_res, lat_res, round_digs, foot_dim_lonxlat,
       
       #flux start
       flux_year_start, flux_month_start, flux_day_start,
       flux_hour_start, flux_min_start,
       
       #flux end
       flux_year_end, flux_month_end, flux_day_end,
       flux_hour_end, flux_min_end,
       
       #flux res
       flux_t_res,
       
       # obs/receptor start
       obs_year_start, obs_month_start, obs_day_start,
       obs_hour_start, obs_min_start,
       
       # obs/receptor end
       obs_year_end, obs_month_end, obs_day_end,
       obs_hour_end, obs_min_end,
       
       #obs resolution
       obs_t_res,
       
       subset_hours_utc,
       aggregate_obs,
       min_agg_obs,
       ntimes,
       tstart,
       tstop,
       flux_units,
       ls, lt,
       
       #save here
       file = file.path('X_config.RData')
  )
  
} #closes function