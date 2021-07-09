prepare.parallelization <- function(homedir, input.extension, workdir, footprint.dirs,
                                    odiac.path, vulcan.path, is.annual, vulcan.sector.path, smurf.path,
                                    use.year, carma.path, edgar.path, edgar.downscaling, edgar.sector.path,
                                    output.directory, site, xmin, xmax, ymin, ymax, emissions.category,
                                    api.key, included.errors, prior.emissions, truth.emissions,
                                    background.emissions, prior.uncertainty, bio.emissions, prior.lonlat,
                                    background.lonlat, observation.values, background.values,
                                    ls = ls, lt = lt, obs.error, lon_res = lon_res, lat_res = lat_res,
                                    flux_units = flux_units) {
  
  #' Get a list of the X-STILT *footprint* directories.
  #' The basename of these directories will become the
  #' output directories for the Bayesian inversion code.
  #' First, the appropriate directories need to be constructed
  #' and the Bayesian Inversion scripts added. This is done using
  #' the `output.directory` variable. In each directory, three more
  #' directories must be included:
  #' *footprints*, *include*, *inversion*
  out.dirs <- file.path(output.directory, basename(footprint.dirs))
  necessary.dirs <- file.path(out.dirs,
                              rep(c('footprints', 'include', 'inversion'),
                                  each = length(out.dirs)))
  
  #create the directories
  for(i in 1:length(necessary.dirs)) {
    dir.create(necessary.dirs[i], recursive = T, showWarnings = F)
  }
  
  #add the Bayesian Inversion scripts to the new directories
  #keeping it simple... move the *inversion* scripts first
  for(i in 1:length(out.dirs)) {
    file.copy(from = list.files('r/src/inversion_scripts/inversion',
                                full.names = T, pattern = '*.r$'),
              to = file.path(out.dirs[i], 'inversion'),
              overwrite = TRUE)
  }
  
  #now, move the *src* scripts into the inversion directories
  for(i in 1:length(out.dirs)) {
    
    #check for the existence of the src directory
    src.filepath <- file.path(out.dirs[i], 'inversion', 'src')
    if(!dir.exists(src.filepath)) dir.create(src.filepath)
    
    file.copy(from = list.files('r/src/inversion_scripts/inversion/src',
                                full.names = TRUE),
              to = src.filepath, overwrite = TRUE)
    
  }
  
  #finally, move the included.errors.csv file into the directory for R_config
  for(i in 1:length(out.dirs)) {
    inversion.filepath <- file.path(out.dirs[i], 'inversion')
    file.copy(from = 'ext/included.errors.csv', to = inversion.filepath,
              overwrite = TRUE)
  }
  
  #construct the dataframe for BayesaianInversion.r here
  args.out <- data.frame(site = site,
                         workdir = workdir,
                         api.key = api.key,
                         included.errors = included.errors,
                         footprint.dirs = footprint.dirs,
                         odiac.path = odiac.path,
                         vulcan.path = vulcan.path,
                         is.annual = is.annual,
                         vulcan.sector.path = vulcan.sector.path,
                         edgar.path = edgar.path,
                         edgar.downscaling = edgar.downscaling,
                         carma.path = carma.path,
                         smurf.path = smurf.path,
                         use.year = use.year,
                         output.directory = output.directory,
                         emissions.category = emissions.category,
                         edgar.sector.path = edgar.sector.path,
                         xmin = xmin, xmax = xmax,
                         ymin = ymin, ymax = ymax,
                         prior.emissions = prior.emissions,
                         truth.emissions = truth.emissions,
                         prior.lonlat = prior.lonlat,
                         background.emissions = background.emissions,
                         background.lonlat = background.lonlat,
                         bio.emissions = bio.emissions,
                         prior.uncertainty = prior.uncertainty,
                         observation.values = observation.values,
                         background.values = background.values,
                         lon_res = lon_res, lat_res = lat_res,
                         ls = ls, lt = lt, obs.error = obs.error,
                         flux_units = flux_units)
  return(args.out)
  
}
