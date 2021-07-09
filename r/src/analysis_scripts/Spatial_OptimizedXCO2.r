if(!file.exists('Out/Spatial_XCO2diff.csv')) {

  #Spatial_XCO2_Optimization
  SpatialData <- data.frame(matrix(NA, nrow = 0, ncol = 5))
  names(SpatialData) <- c('time', 'lon', 'lat',
                          'prior.xco2', 'posterior.xco2')
    
  for(i in 1:length(dir.list)) {
    
    ###############################################################
    #' Construct the *include*, *out*, and *footprint* paths
    #' corresponding to the inversion's output directories.
    #' Finally, read in the `prior` and `posterior` data.
    ###############################################################
    include.path <- file.path(dir.list[i], 'include')
    out.path <- file.path(dir.list[i], 'inversion/out')
    foot.path <- file.path(dir.list[i], 'footprints')
    
    #list all of the *.nc footprint files
    footprints <- list.files(foot.path, recursive = TRUE,
                             full.name = TRUE, pattern = '*.nc')
    
    #get the prior and posterior *.nc files
    prior_emiss <- brick(file.path(include.path, 'prior_emiss.nc'))
    posterior_emiss <- brick(file.path(out.path, 'posterior.nc'))
    ###############################################################
    
    
    
    ###############################################################
    #' Loop through each footprint and multiply by the prior and
    #' posterior rasters.
    ###############################################################
    for(j in 1:length(footprints)) {
      
      #read in the footprint raster.
      footprint.raster <- brick(footprints[j])
      
      #Get the date, time, lon, and lat from the file name
      tmp.data <- 
        str_split_fixed(gsub('.nc', '', basename(footprints[j])),
                        pattern = '_', n = 3)
      add.line <- data.frame(time = tmp.data[1,1],
                             lon = as.numeric(tmp.data[1,2]),
                             lat = as.numeric(tmp.data[1,3]),
                             prior.xco2 = NA,
                             posterior.xco2 = NA)
      
      prior.xco2 <- posterior.xco2 <- NULL
      for(k in 1:nlayers(footprint.raster)) {
        
        # get the footprint layer
        footprint.layer <- footprint.raster[[k]]
        
        #' Get the name of the footprint file.
        #' Use this name to select the layers from
        #' the other raster files. Each file's layers
        #' may be in a different order than the 
        #' footprint's.
        layer.name <- names(footprint.layer)
        
        # get the prior layer
        eval(parse(text = paste0('prior.layer <- prior_emiss$',
                                 layer.name)))
        
        # get the posterior layer
        eval(parse(text = paste0('posterior.layer <- posterior_emiss$',
                                 layer.name)))
        
        # add the hourly xco2 contribution to the appropriate vector
        # The footprint layer will likely be larger than the prior layer
        prior.xco2[k] <- 
          suppressWarnings(cellStats(footprint.layer*prior.layer, sum))
        posterior.xco2[k] <-
          suppressWarnings(cellStats(footprint.layer*posterior.layer, sum))
        
      } # close the hourly xco2 list
      
      # sum up the hourly constributions
      add.line$prior.xco2 <- sum(prior.xco2)
      add.line$posterior.xco2 <- sum(posterior.xco2)
      
      # add the current calculations to the complete dataframe
      SpatialData <- rbind(SpatialData, add.line)
      ###############################################################
      
    } # close the footprint list
  } # close the directory list
  
  
  ###############################################################
  #' Since it takes so long to make the `SpatialData` object,
  #' a *.csv copy is saved to the output directory.
  ###############################################################
  # convert the time to POSIXt and change tz to local
  SpatialData$time <- as.POSIXct(SpatialData$time,
                                 format = '%Y%m%d%H%M',
                                 
                                 tz = 'UTC')
  attr(SpatialData$time, 'tzone') <- local.tz
  
  write.csv(SpatialData, file = 'Out/Spatial_XCO2diff.csv',
            row.names = FALSE)
  ###############################################################
  
  # closes the conditional (Read Spatial_XCO2diff.csv if it exists)
} else {SpatialData <- read.csv('Out/Spatial_XCO2diff.csv')} 

# Save as a character string
SpatialData$time <- as.character(SpatialData$time)

###############################################################
#' Now, create some additional, non-spatial plots that may be 
#' beneficial in the analysis process.
###############################################################
#determine the mean XCO2 in the prior and posterior data
mean.SpatialData <- aggregate(SpatialData[,4:5],
                              list(round(hour(SpatialData$time) +
                                           minute(SpatialData$time)/60,
                                   2)),
                              mean)
mean.SpatialData$Group.1 <- as.character(mean.SpatialData$Group.1)

#create a new SpatialData dataframe for plotting purposes
tmp1.SpatialData <- data.frame(SpatialData[,1:4], 'prior')
names(tmp1.SpatialData) <-
  c('time', 'lon', 'lat', 'xco2', 'Emission Inventory')
tmp2.SpatialData <- data.frame(SpatialData[,c(1:3,5)], 'posterior')
names(tmp2.SpatialData) <- names(tmp1.SpatialData)

#combine both dataframe
tmp.SpatialData <- rbind(tmp1.SpatialData, tmp2.SpatialData)
tmp.SpatialData$hour.string <- as.character(
  round(hour(tmp.SpatialData$time) + minute(tmp.SpatialData$time)/60, 2)
)

#reorder the hour.strings as factors
tmp.SpatialData$hour.string <-
  factor(tmp.SpatialData$hour.string, levels = mean.SpatialData$Group.1)


Boxplots_XCO2 <-
  ggplot() +
  ggtitle(expression(paste('Hourly Changes in ', Delta, 'XCO'[2]))) +
  labs(subtitle = paste0('Site: ', site),
       caption = paste0('Black Line: Prior Mean', '\n',
                        'Red Line: Posterior Mean')) +
  xlab('Hour of Day [hr]') +
  ylab(expression(paste(Delta, 'XCO'[2]))) +
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 'dashed') +
  geom_boxplot(data = tmp.SpatialData,
               aes(x = hour.string,
                   y = xco2, fill = `Emission Inventory`)) +
  stat_summary(data = subset(tmp.SpatialData,
                             `Emission Inventory` == 'prior'),
               aes(x = hour.string,
                   y = xco2, group = 1), fun = mean,
               geom = 'line', linetype = 'dashed', color = 'black') +
  stat_summary(data = subset(tmp.SpatialData,
                             `Emission Inventory` == 'posterior'),
               aes(x = hour.string,
                   y = xco2, group = 1), fun = mean,
               geom = 'line', linetype = 'dashed', color = 'red') +
  theme_classic() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')

# save the plot
ggsave(Boxplots_XCO2,
       filename = file.path('Out', 'Hourly_OptimizedXCO2.jpg'),
       device = 'jpg', width = 8, height = 6, units = 'in')

###############################################################

# plot the data
Optimized_XCO2.plot <- 
  ggmap(map.10) +
  labs(subtitle = '[Local Time]') +
  ggtitle(expression(paste('Optimized XCO'[2], ' Over the L.A. Basin'))) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_point(data = SpatialData,
             aes(x = lon, y = lat, fill = posterior.xco2 - prior.xco2),
             shape = 23) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                       name = expression(paste(Delta, 'XCO'[2]))) +
  facet_wrap(. ~ strftime(time, tz = local.tz), ncol = 3) +
  theme(text = element_text(size = 9),
       plot.title = element_text(hjust = 0.5),
       legend.position = 'bottom',
       legend.key.width = unit(0.7, 'in'))
  
# save the plot
ggsave(Optimized_XCO2.plot,
       filename = file.path('Out', 'Spatial_XCO2diff.jpg'),
       device = 'jpg', width = 8, height = 9, units = 'in')
