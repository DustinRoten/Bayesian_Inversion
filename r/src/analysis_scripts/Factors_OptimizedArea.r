library(ggplot2); library(patchwork)
library(raster); library(geodist)
library(lubridate); library(plot3D)

#determine the center of the domain of interest
mean.x <- (xmin + xmax)/2
mean.y <- (ymin + ymax)/2

if(!exists('Error.Reduction')) {
  
  Error.Reduction <-
    data.frame(matrix(NA, nrow = length(dir.list), ncol = 9))
  names(Error.Reduction) <- c('Soundings', 'mean.Time', 'mean.Hour',
                              'mean_km_covered', 'mean_xco2', 'sd_xco2',
                              'mean_DiffFlux', 'sd_DiffFlux',
                              'Perc_Uncert_Red')
  for(i in 1:length(dir.list)) {
    
    ###############################################################
    #' Construct the *include* and *out* paths corresponding
    #' to the inversion's output directories
    ###############################################################
    include.path <- file.path(dir.list[i], 'include')
    out.path <- file.path(dir.list[i], 'inversion/out')
    ###############################################################
    
    
    
    ###############################################################
    #' For the first calculation, find the difference between
    #' the `prior` and `posterior emissions`. Include the areas
    #' of cells where the difference is greater than 1umol/m2/s.
    ###############################################################
    prior_emiss <- brick(file.path(include.path, 'prior_emiss.nc'))
    posterior_emiss <- brick(file.path(out.path, 'posterior.nc'))
    obs <- readRDS(file.path(include.path, 'obs.rds'))
    area.km2 <- NULL; diff.flux <- NULL
    for(j in 1:nlayers(prior_emiss)) {
      #obtain a single prior and posterior layer
      prior_layer <- prior_emiss[[j]]
      eval(parse(text = paste0('posterior_layer <- posterior_emiss$',
                               names(prior_layer))))
      
      #calculate the difference between the layers
      diffs <- posterior_layer - prior_layer
      bin.diffs <- diffs
      bin.diffs[abs(bin.diffs) >= 1] <- 1
      area.bin.diffs <- raster::area(bin.diffs)
      
      #Add values to vectors here
      area.km2[j] <- cellStats(area.bin.diffs*bin.diffs, sum)
      diff.flux <- c(diff.flux, values(diffs))
      
    }
    
    max.km2 <- max(area.km2); mean.diff <- mean(diff.flux)
    sd.diff <- sd(diff.flux)
    
    #add the mean area optimized
    Error.Reduction$mean_km_covered[i] <- max.km2
    Error.Reduction$mean_DiffFlux[i] <- mean.diff
    Error.Reduction$sd_DiffFlux[i] <- sd.diff
    ###############################################################
    
    
    
    ###############################################################
    #' Read in the Percent Uncertainty Reduction
    ###############################################################
    perc_unc_red <- readRDS(file.path(out.path, 'perc_unc_red.rds'))
    footprint.dirs <- list.files(file.path(dir.list[i], 'footprints'),
                                 full.names = TRUE)
    num_soundings <- length(footprint.dirs)
    
    Error.Reduction$Soundings[i] <- num_soundings
    Error.Reduction$Perc_Uncert_Red[i] <- perc_unc_red
    ###############################################################
    
    
    
    ###############################################################
    #' Determine the mean time of the SAM observation.
    #' Lon/lat values are included here for future uses.
    #' to each sounding that makes up a SAM.
    ###############################################################
    lon.lat.time <- data.frame(matrix(NA, nrow = 0, ncol = 3))
    names(lon.lat.time) <- c('time', 'lon', 'lat')
    for(j in 1:length(footprint.dirs)) {
      #get the spatial information from the file names
      file.name <- list.files(footprint.dirs[j])
      add.line <- str_split_fixed(file.name, pattern = '_', n = 3)
      #create a temporary line to add to the dataframe
      tmp_df <- data.frame(time = as.POSIXct(add.line[1,1],
                                             format = '%Y%m%d%H%M',
                                             tz = 'UTC'),
                           lon = as.numeric(add.line[1,2]),
                           lat = as.numeric(gsub('.nc', '', add.line[1,3])))
      #add to the dataframe
      lon.lat.time <- rbind(lon.lat.time, tmp_df)
    }
    ###############################################################
    
    #convert to local time
    tmp.time <- mean(lon.lat.time$time)
    attr(tmp.time, 'tzone') <- local.tz
    
    Error.Reduction$mean.Time[i] <- tmp.time
    Error.Reduction$mean.Hour[i] <- hour(tmp.time) + minute(tmp.time)/60
    
    #include mean and sd xco2
    Error.Reduction$mean_xco2[i] <- mean(as.numeric(obs[,3]))
    Error.Reduction$sd_xco2[i] <- sd(as.numeric(obs[,3]))
    
  }
} #closes if statement (Error.Reduction exists)

Error.Reduction$mean.Time <- as.POSIXct(Error.Reduction$mean.Time,
                                        origin = '1970-01-01',
                                        tz = local.tz)

##### Make Plots! #####

###################################
### 2 Plots - Soundings vs Time ###
###################################
### Time of Year
linear.model.1 <- lm(Soundings ~ mean.Time, data = Error.Reduction)
Soundings_vs_ToY <-
  ggplot() +
  ggtitle('Number of Soundings vs. Time of Year') +
  labs(subtitle = paste0('Site: ', site, '\n', '[Running Average]')) +
  xlab('Time of Year') +
  ylab('Number of Soundings') +
  geom_smooth(data = Error.Reduction,
              aes(x = mean.Time, y = Soundings),
              color = 'red', size = 0.5,
              linetype = 'dashed', se = FALSE) +
  geom_point(data = Error.Reduction,
             aes(x = mean.Time, y = Soundings)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        axis.text = element_text(size = 10))

### Time of Day
linear.model.2 <- lm(Soundings ~ mean.Hour, data = Error.Reduction)
Soundings_vs_ToD <-
  ggplot() +
  ggtitle('Number of Soundings vs. Time of Day') +
  labs(subtitle = paste0('Site: ', site, '\n', '[Ordinary Least Squares Regression]')) +
  xlab('Time of Day [hr]') +
  ylab('Number of Soundings') +
  geom_smooth(data = Error.Reduction,
              aes(x = mean.Hour, y = Soundings),
              color = 'red', size = 0.5,
              linetype = 'dashed', method = 'lm') +
  geom_point(data = Error.Reduction,
             aes(x = mean.Hour, y = Soundings)) +
  geom_text(aes(x = Inf, y = Inf, vjust = 1, hjust = 1,
                label = paste0(' Slope = ',
                               round(linear.model.2$coefficients[2], 3), '\n',
                               ' Y-Intercept = ',
                               round(linear.model.2$coefficients[1], 3), '\n',
                               ' R-Squared = ',
                               round(summary(linear.model.2)$r.squared, 3), '\n',
                               ' P-Value = ',
                               round(summary(linear.model.2)$coefficients[2,4], 3)
                ))) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        axis.text = element_text(size = 10))
###################################



##########################
### Area vs. Soundings ###
##########################
### Next, Area vs. SAM Size
Linearized.Data <-
  data.frame(mod.soundings = 1/Error.Reduction$Soundings,
             mod.area = 1/Error.Reduction$mean_km_covered)

linear.model.3 <-
  lm(mod.area ~ mod.soundings, data = Linearized.Data)

#calculate alpha
alpha.value <- 1/linear.model.3$coefficients[1]
names(alpha.value) <- 'alpha.value'

#calculate beta
beta.over.alpha <- linear.model.3$coefficients[2]
beta.value <- alpha.value*beta.over.alpha

#expression for ggplot label
my_exp <-
  expression(paste('Saturation Growth Rate Equation: y=',
                   alpha, frac(N,beta + N)))

Area_vs_SAM <- 
  ggplot() +
  ggtitle('Mean Area Optimized vs. SAM Size') +
  labs(subtitle = paste0('Site: ', site, '\n', '[Linearized SGR Equation]')) +
  xlab(expression(paste(frac(1,N), '   [Soundings'^-1, ']'))) +
  ylab(expression(paste(frac(1, Area), '   [km'^-2, ']'))) +
  geom_smooth(data = Linearized.Data,
              aes(x = mod.soundings, y = mod.area),
              method = 'lm', color = 'red', size = 0.5,
              linetype = 'dashed') +
  geom_point(data = Linearized.Data,
             aes(x = mod.soundings, y = mod.area)) +
  geom_text(aes(x = -Inf, y = Inf, vjust = 1, hjust = 0,
                label = paste0(' Slope = ',
                               round(linear.model.3$coefficients[2], 3), '\n',
                               ' Y-Intercept = ',
                               round(linear.model.3$coefficients[1], 3), '\n',
                               ' R-Squared = ',
                               round(summary(linear.model.3)$r.squared, 3), '\n',
                               ' P-Value = ',
                               round(summary(linear.model.3)$coefficients[2,4], 3), '\n',
                               ' Alpha = ',
                               round(alpha.value, 1), '\n',
                               ' Beta = ',
                               round(beta.value, 1)
                ))) + 
  geom_text(aes(x = Inf, y = -Inf, vjust = 0, hjust = 1),
            label = my_exp) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        axis.text = element_text(size = 10))

### Area covered vs Time of Day
linear.model.4 <- lm(mean_km_covered ~ mean.Hour, data = Error.Reduction)
Area_vs_ToD <-
  ggplot() +
  ggtitle('Max Area Optimized vs. Time of Day') +
  labs(subtitle = paste0('Site: ', site, '\n', '[Ordinary Least Squares Regression]')) +
  xlab('Time of Day [hr]') +
  ylab(expression(paste('Max Area Optimized [km'^2, ']'))) +
  geom_smooth(data = Error.Reduction,
              aes(x = mean.Hour, y = mean_km_covered),
              color = 'red', size = 0.5,
              linetype = 'dashed', method = 'lm') +
  geom_point(data = Error.Reduction,
             aes(x = mean.Hour, y = mean_km_covered)) +
  geom_text(aes(x = Inf, y = Inf, vjust = 1, hjust = 1,
                label = paste0(' Slope = ',
                               round(linear.model.4$coefficients[2], 3), '\n',
                               ' Y-Intercept = ',
                               round(linear.model.4$coefficients[1], 3), '\n',
                               ' R-Squared = ',
                               round(summary(linear.model.4)$r.squared, 3), '\n',
                               ' P-Value = ',
                               round(summary(linear.model.4)$coefficients[2,4], 3)
                ))) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 10),
        axis.text = element_text(size = 10))

### Make the Quad Plot
quad.plot <-
  (Soundings_vs_ToY | Soundings_vs_ToD) / (Area_vs_SAM | Area_vs_ToD) +
  plot_annotation(title = 'Factors Influencing Area of Domain Optimization',
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 15)),
                  tag_levels = 'A')

# save the plot
ggsave(quad.plot, filename = file.path('Out', 'Factors_OptimizedArea.jpg'),
       device = 'jpg', width = 10, height = 9, units = 'in')

