#inversion.analysis.3.r
#Run Bayesian inversion
#grab all of the csv files
TEST <- FALSE
options(stringsAsFactors = FALSE)
#################################
##### Setup user parameters #####
#################################
homedir <- '/uufs/chpc.utah.edu/common/home'
input.extension <- 'lin-group11/Roten_InputData'
workdir <- '~/Sectoral_Analysis'
setwd(workdir); source('r/dependencies.r')
library(patchwork); library(geodist); library(RNetCDF)

#site name
site <- 'Los Angeles'

#footprints data path
footprint.type <- 'ODIAC'
footprint.dirs <- list.files(file.path(homedir,
                                       'lin-group11/XCO2_Climatology',
                                       footprint.type),
                             full.names = TRUE, pattern = 'out_')

#ODIAC data path
odiac.path <-
  file.path(homedir, input.extension, 'ODIAC', 'ODIAC2020', '2019')

#Vulcan data path (hourly or annual emissions?)
#sector dataframe (Vulcan)
vulcan.path <- file.path(homedir, input.extension,
                         'Vulcan3.0/hourly')
is.annual <- FALSE
vulcan.sector.path <- 'ext/defined_vulcan_sectors.csv'

#SMUrF data path
smurf.path <- file.path(homedir, input.extension, 'SMUrF')
use.year <- 2019 #not all of 2020 is available

#CARMA dataset path
carma.path <- 'ext/CARMA/Plant.csv'

#EDGAR data path
edgar.path <- file.path(homedir, 'lin-group11/Roten_InputData',
                        'EDGARv5', '2018')

#EDGAR temporal downscaling
edgar.downscaling <- 'ext/EDGAR_TemporalProfiles'

#sector dataframe (EDGAR)
edgar.sector.path <- 'ext/defined_edgar_sectors.csv'

#identify the location of the output
output.directory <- file.path(homedir, 'lin-group11/Bayesian_Inversion')

#domain of interest
site <- 'Los Angeles'
xmin <- -118.6878; xmax <- -117.0617
ymin <- 33.58086; ymax <- 34.34599

#resolution
lon_res <- 1/120; lat_res <- 1/120

#scaling parameters
ls <- 6; lt <- 1

#observation error (in ppm)
obs.error <- 0.23

#flux units
flux_units <- 'umol/(m2 s)'

#identify the sector of interest here
#use NA for the total XCO2 amount
#use 'Biosphere' for biospheric flux
#must be in both EDGAR and Vulcan sectors!
emissions.category <- NA

#register the api key for google maps to work
api.key <-
  read.csv('insert_ggAPI.csv', header = FALSE, col.names = 'api.key')
register_google(key = api.key)

#simulated XCO2 observation errors
#construct a dataframe with the name of each error source
#second column is the mean value of the error (or NA if it is a %)
#third column is the standard deviation of the eorror (or NA if it is a %)
#fourth column is the percent error (other numeric columns must be NA)
included.errors <- 'ext/included.errors.csv'

#set up file names for required NetCDF files
prior.emissions <- 'prior_emiss.nc'
truth.emissions <- 'truth_emiss.nc'
background.emissions <- 'outer_emiss.nc'
prior.uncertainty <- 'prior_uncert.nc'
bio.emissions <- 'bio_flux.nc'

#set up file names for required RDS files
prior.lonlat <- 'lonlat_domain.rds'
background.lonlat <- 'lonlat_outer.rds'
observation.values <- 'obs.rds'
background.values <- 'background.rds'

#include SLURM options
slurm_options <- list(time = '96:00:00',
                      account = 'lin-np',
                      partition = 'lin-np')
#################################
#################################
#################################

#create a dataframe to be sent to SLURM
p.table <-
  prepare.parallelization(homedir = homedir,
                          input.extension = input.extension,
                          workdir = workdir,
                          api.key = api.key,
                          footprint.dirs = footprint.dirs,
                          odiac.path = odiac.path,
                          vulcan.path = vulcan.path,
                          is.annual = is.annual,
                          vulcan.sector.path = vulcan.sector.path,
                          smurf.path = smurf.path,
                          use.year = use.year,
                          carma.path = carma.path,
                          edgar.path = edgar.path,
                          edgar.downscaling = edgar.downscaling,
                          edgar.sector.path = edgar.sector.path,
                          output.directory = output.directory,
                          site = site,
                          xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          emissions.category = emissions.category,
                          included.errors = included.errors,
                          prior.emissions = prior.emissions,
                          truth.emissions = truth.emissions,
                          background.emissions = background.emissions,
                          prior.uncertainty = prior.uncertainty,
                          bio.emissions = bio.emissions,
                          prior.lonlat = prior.lonlat,
                          background.lonlat = background.lonlat,
                          observation.values = observation.values,
                          background.values = background.values,
                          lon_res = lon_res, lat_res = lat_res,
                          ls = ls, lt = lt, obs.error,
                          flux_units = flux_units)

#for testing
if(TEST == TRUE) p.table <- p.table[1,]

####################################################
##### Begin the Inversion Process for each SAM #####
####################################################
#run the Bayesian Inversion in parallel!
if(nrow(p.table) > 1) {
  
  # run jobs in parallel with SLURM
  rslurm::slurm_apply(f = BayesianInversion,
                      params = p.table,
                      jobname = 'Inversion',
                      nodes = 1,
                      cpus_per_node = nrow(p.table),
                      pkgs = 'base',
                      slurm_options = slurm_options)
  
} else if (nrow(p.table == 1)) {
  
  # run one job on an interactive node
  BayesianInversion(site = p.table$site,
                    workdir = p.table$workdir,
                    api.key = p.table$api.key,
                    included.errors = p.table$included.errors,
                    footprint.dirs = p.table$footprint.dirs,
                    odiac.path = p.table$odiac.path,
                    vulcan.path = p.table$vulcan.path,
                    is.annual = p.table$is.annual,
                    vulcan.sector.path = p.table$vulcan.sector.path,
                    edgar.path = p.table$edgar.path,
                    edgar.downscaling = p.table$edgar.downscaling,
                    carma.path = p.table$carma.path,
                    smurf.path = p.table$smurf.path,
                    use.year = p.table$use.year,
                    output.directory = p.table$output.directory,
                    emissions.category = p.table$emissions.category,
                    edgar.sector.path = p.table$edgar.sector.path,
                    xmin = p.table$xmin,
                    xmax = p.table$xmax,
                    ymin = p.table$ymin,
                    ymax = p.table$ymax,
                    prior.emissions = p.table$prior.emissions,
                    prior.lonlat = p.table$prior.lonlat,
                    background.emissions = p.table$background.emissions,
                    background.lonlat = p.table$background.lonlat,
                    bio.emissions = p.table$bio.emissions,
                    prior.uncertainty = p.table$prior.uncertainty,
                    observation.values = p.table$observation.values,
                    background.values = p.table$background.values,
                    lon_res = p.table$lon_res,
                    lat_res = p.table$lat_res,
                    ls = p.table$ls,
                    lt = p.table$lt,
                    obs.error = p.table$obs.error,
                    flux_units = p.table$flux_units)
}