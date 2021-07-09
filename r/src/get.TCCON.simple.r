#' D. Roten (01/2021)
#' This "simple" version of `get.TCCON()` requires the 
#' exact file path to the TCCON data (*.csv format). Using
#' this file and a provided range (+/- hrs), an average 
#' is calculated.
get.TCCON.simple <- function(tccon.file = NULL,
                             date.time = NULL, range = NULL) {
  
  TCCON.data <- read.csv(tccon.file)
  TCCON.data$times <- as.POSIXlt(TCCON.data$times,
                                 origin = '1970-01-01',
                                 tz = 'UTC')
  
  sub.TCCON.data <- subset(TCCON.data,
                           (times > (date.time-range*3600)) &
                             (times < (date.time+range*3600)))
  
  TCCON.xco2 <- mean(sub.TCCON.data$xco2)
  TCCON.xco2.uncert <-
    sd(sub.TCCON.data$xco2)/sqrt(nrow(sub.TCCON.data))
  
  TCCON.df <- data.frame(cbind(TCCON.xco2,
                               TCCON.xco2.uncert))
  
  return(TCCON.df)
  
}