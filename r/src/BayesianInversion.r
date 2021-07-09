BayesianInversion <- function(site, workdir, api.key, included.errors, footprint.dirs, odiac.path,
                              vulcan.path, is.annual, vulcan.sector.path, edgar.path, edgar.downscaling,
                              carma.path, smurf.path, use.year, output.directory, emissions.category,
                              edgar.sector.path, xmin, xmax, ymin, ymax, prior.emissions, truth.emissions,
                              prior.lonlat, background.emissions, background.lonlat, bio.emissions,
                              prior.uncertainty, observation.values, background.values,
                              lon_res, lat_res, ls, lt, obs.error, flux_units) {
  
  setwd(workdir); source('r/dependencies.r')
  library(patchwork); library(geodist); library(RNetCDF)
  library(ggmap)
  
  #if needed, add api.key after workdir in the arguments
  register_google(key = api.key)
  
  #read in external data (from /ext)
  included.errors <- read.csv(included.errors, header = TRUE)
  edgar.sector_df <- read.csv(edgar.sector.path)
  vulcan.sector_df <- read.csv(vulcan.sector.path)
  
  #generate prior extent (inner domain)
  nc.extent <- extent(xmin, xmax, ymin, ymax)
  
  for(i in 1:length(footprint.dirs)) {
    
    #' First, the appropriate directories need to be constructed
    #' and the Bayesian Inversion scripts added. This is done using
    #' the `output.directory` variable. In each directory, three more
    #' directories must be included:
    #' *footprints*, *include*, *inversion*
    out.dirs <- file.path(output.directory, basename(footprint.dirs[i]))
    necessary.dirs <- file.path(out.dirs,
                                rep(c('footprints', 'include', 'inversion'),
                                    each = length(out.dirs)))

    #create the directories
    for(dir in 1:length(necessary.dirs))
      dir.create(necessary.dirs[dir], recursive = T, showWarnings = F)
    
    #create the operator.dir (footprints)
    operator.dir <- file.path(out.dirs, 'footprints')
    obs.dir <- file.path(out.dirs, 'include')
    
    #' Compile the list of footprints here. The standard output of
    #' X-STILT includes a `footprints` folder of symbolic links to
    #' the actual footprint files. If this directory is present, 
    #' the file path will be automatically updated. If it is NOT 
    #' present, the script will search for footprints within the 
    #' original specified directory.
    listed.dir <- list.files(footprint.dirs[i], full.names = TRUE)
    idx <- which(basename(listed.dir) == 'footprints')

    #update directory as needed
    if(length(idx) != 0)
      footprint.list <- list.files(listed.dir[idx],
                                   full.names = TRUE, pattern = '.nc')
    if(length(idx) == 0)
      footprint.list <- list.files(footprint.dirs[i],
                                   full.names = TRUE, pattern = '.nc')
    
    ######################################################
    #' The maximum extent required for the background
    #' emission inventory file must be determined. The
    #' low-tech approach here simply reads in each 
    #' footprint file and records the maximum lons/lats.
    #' Additionally, averaged timesteps will be generated.
    ######################################################
    min.lon <- NULL; max.lon <- NULL
    min.lat <- NULL; max.lat <- NULL
    start.times <- NULL; timesteps <- NULL
    for(j in 1:length(footprint.list)) {
      
      #read in the raster, remove all zero weights
      footprint <- brick(footprint.list[j])
      footprint_df <- raster::as.data.frame(footprint, xy = TRUE)
      summed.flux.values <- rowSums(footprint_df[,3:ncol(footprint_df)])
      tmp_df <- data.frame(footprint_df[,1:2], summed.flux.values)
      tmp_df <- subset(tmp_df, summed.flux.values > 0)
      
      #record the number of layers
      timesteps[j] <- nlayers(footprint)
      
      #Determine the min/max values here
      min.lon[j] <- min(tmp_df$x); max.lon[j] <- max(tmp_df$x)
      min.lat[j] <- min(tmp_df$y); max.lat[j] <- max(tmp_df$y)
      
      #' Record the start time of each footprint.
      #' THE START TIME MUST BE IN THE NAME OF THE FILE!
      start.times[j] <-
        unlist(strsplit(basename(footprint.list[j]), split = '_'))[1]
      
      cat(paste0('Constraining footprint domain: ',
                 round(j/length(footprint.list), 2), '%'), '\r')
      
    }; remove('footprint'); remove('tmp_df')
    
    #Overwrite the unnecessary vectors and replace with a single
    #min/max value.
    #' Testing Note: For SAM_2020022419, the values are as follows:
    #' min.lon <- -121.6875; max.lon <- -115.4708
    #' min.lat <- 32.79583; max.lat <- 41.10417
    min.lon <- min(min.lon); max.lon <- max(max.lon)
    min.lat <- min(min.lat); max.lat <- max(max.lat)
    outer.extent <- extent(min.lon, max.lon, min.lat, max.lat)
    max.layers <- max(timesteps)
    
    #construct the averaged timesteps here
    #HOURLY steps are assumed!
    #strip off seconds!!
    pos.start.times <- as.POSIXlt(start.times,
                                  format = '%Y%m%d%H%M',
                                  tz = 'UTC')
    avg.start.time <- mean(pos.start.times)
    avg.start.time <- avg.start.time - second(avg.start.time)
    
    #format the average start time into a YYYYMMDDHHMM timestamp
    timestr <- strftime(avg.start.time,
                        format = '%Y%m%d%H%M',
                        tz = 'UTC')
    
    #build the sequence here
    step.in.seconds <- 3600
    avg.timesteps <-
      seq(as.numeric(avg.start.time),
          as.numeric(avg.start.time) - step.in.seconds*max.layers,
          -1*step.in.seconds)
    
    #remove the start time
    #(its not included in the footprint timesteps)
    timestep.list <- as.POSIXlt(avg.timesteps,
                                origin = '1970-01-01',
                                tz = 'UTC')[2:length(avg.timesteps)]
    
    #########################################
    #' Construct the modified footprint files
    #########################################
    #' When a footprint is generated, the timesteps will be
    #' mapped to the averaged timesteps from above. The footprint
    #' will also be cropped according to the outer domain that
    #' was constrained above.
    make_foots_ncdf4(footprint.list = footprint.list,
                     timestr = timestr,
                     averaged.timesteps = timestep.list,
                     operator.dir = operator.dir,
                     nc.extent = outer.extent)
    gc() #garbage collection
    #########################################
    
    
    ########################################################
    #' Make the `prior_emiss.nc` and `outer_emiss.nc` files.
    ########################################################
    make_prior_ncdf4(site, odiac.path, edgar.path,
                     category = emissions.category,
                     sector.list = edgar.sector_df,
                     times = timestep.list,
                     inner.extent = nc.extent,
                     full.extent = outer.extent,
                     downscaling.extension = edgar.downscaling,
                     carma.path, obs.dir,
                     inner.name = prior.emissions,
                     outer.name = background.emissions)
    gc() #garbage collection
    ########################################################
    
    
    ####################################################
    #' Make the `bio_flux.nc` file using `SMUrF` output.
    ####################################################
    #Grab the outer emissions file to use as a reference.
    #Any layer will do. No need to read the entire file.
    outer.domain <- brick(file.path(obs.dir, background.emissions))[[1]]
    
    #make the layered (NetCDF) file for the biospheric flux (from SMUrF)
    make_bio_ncdf4(SMUrF.path = smurf.path,
                   times = timestep.list,
                   match.year = use.year,
                   match.raster = outer.domain,
                   bio.output = file.path(obs.dir,
                                          bio.emissions))
    gc()
    ####################################################
    
    
    #############################################################
    #' Next, the files `lonlat_domain.rds` and `lonlat_outer.rds`
    #' will be generated by using the output from the previous
    #' function: `make_prior_ncdf4()`.
    #############################################################
    #get the values for the inner domain
    for.lonlat_domain <- nc_open(file.path(obs.dir, prior.emissions))
    Lon <- ncvar_get(for.lonlat_domain, 'lon')
    Lat <- ncvar_get(for.lonlat_domain, 'lat')
    nc_close(for.lonlat_domain)
    
    #create the matrix of data and save (inner domain)
    lonlat_domain <-
      as.matrix(cbind(rep(Lon, each = length(Lat)), Lat))
    colnames(lonlat_domain) <- c('Lon', 'Lat')
    saveRDS(lonlat_domain,
            file = file.path(obs.dir, prior.lonlat))
    
    #now, the same is done for the outer domain
    for.lonlat_outer <- nc_open(file.path(obs.dir, background.emissions))
    Lon <- ncvar_get(for.lonlat_outer, 'lon')
    Lat <- ncvar_get(for.lonlat_outer, 'lat')
    nc_close(for.lonlat_outer)
    
    #create the matrix of data and save (outer domain)
    lonlat_outer <-
      as.matrix(cbind(rep(Lon, each = length(Lat)), Lat))
    colnames(lonlat_outer) <- c('Var1', 'Var2')
    saveRDS(lonlat_outer,
            file = file.path(obs.dir, background.lonlat))
    #############################################################
    
    
    ################################################################
    #' Construct `emiss_uncert.nc` here by comparing the inner
    #' domain with hourly Vulcan 3.0 estimates (Gurney et al., 2020)
    ################################################################
    make_sigma_ncdf4(prior_emiss = file.path(obs.dir, prior.emissions),
                     vulcan_emiss = vulcan.path, is.annual,
                     sector.list = vulcan.sector_df,
                     emissions.category,
                     sigma.output = file.path(obs.dir, prior.uncertainty))
    gc()
    ################################################################
    
    
    #######################################################
    #' Create the `obs.rds` and `background.rds` files here
    #######################################################
    #' Biospheric emissions can be ignored by setting
    #' `bio.emissions = NA`. Only anthropogenic emissions
    #' will be considered.
    synthetic.obs(footprints = footprint.list,
                  start.time = avg.start.time,
                  errors = included.errors,
                  obs.error = obs.error,
                  inner.emissions = file.path(obs.dir,
                                              truth.emissions),
                  outer.emissions = file.path(obs.dir,
                                              background.emissions),
                  biosphere.emissions = file.path(obs.dir, bio.emissions),
                  obs.output = file.path(obs.dir, observation.values),
                  background.output = file.path(obs.dir,
                                                background.values))
    gc()
    ###############################################
    ###############################################
    ###############################################
    
    
    ########################################################
    #' All of the required files have been generated from
    #' the above code. Now, the inversion process can begin.
    ########################################################
    
    #save to new relative working directory for the inversion
    #change working directory here
    new.workdir <- file.path(out.dirs, 'inversion')
    setwd(new.workdir)
    
    #' set up a new *config.r* file: *X_config.r*
    setup.config(timesteps = timestep.list,
                 src_path = "src/",
                 out_path = "out/",
                 foot_dir = "../footprints/",
                 include_dir = "../include/",
                 clear_H = TRUE,
                 compute_chi_sq = FALSE,
                 include_outer = TRUE,
                 include_bio = TRUE,
                 lonlat_domain_file.NAME = prior.lonlat,
                 lonlat_outer_file.NAME = background.lonlat,
                 prior_file.NAME = prior.emissions, 
                 sigma_file.NAME = prior.uncertainty,
                 outer_file.NAME = background.emissions,
                 bio_file.NAME = bio.emissions,
                 obs_file.NAME = observation.values,
                 bg_file.NAME = background.values,
                 lon_res = lon_res,
                 lat_res = lat_res,
                 foot_dim_lonxlat = T,
                 aggregate_obs = FALSE,
                 min_agg_obs = 2,
                 ls = ls, lt = lt,
                 flux_units = "umol/(m2 s)",
                 save.to = new.workdir)
    
    #start the inversion process
    source('run_all.r')
    
  } #closes footprint.dir loop
} #close function
