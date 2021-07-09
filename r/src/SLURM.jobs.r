SLURM.jobs <- function(partition = NULL) {
  
  if(is.null(partition)) stop('SLURM.jobs: No partition selected...')
  # Raw list of the submissions are included here
  cmd.check <- paste('squeue', '-A', partition, sep = ' ')
  check.subs <- system(cmd.check, intern = TRUE)
  
  # Retrieve the column names from the listed jobs
  split.title <- unlist(strsplit(check.subs[1], split = ''))
  reduced.string <- NULL # Store the column names here
  for(i1 in 1:(length(split.title)-1)) {
    
    if(split.title[i1] != " ") reduced.string[length(reduced.string)+1] <- split.title[i1]
    if(i1 > 1 & split.title[i1] != " " & split.title[i1+1] == " ") reduced.string[length(reduced.string)+1] <- " "
    if(i1 == (length(split.title)-1)) reduced.string[length(reduced.string)+1] <- split.title[i1+1]
    
  }
  
  # Generate data frame
  SLURM_OUT <- data.frame(matrix(NA, nrow = (length(check.subs)-1),
                                 ncol = length(unlist(strsplit(paste(reduced.string, collapse = ''), split = ' ')))))
  names(SLURM_OUT) <- unlist(strsplit(paste(reduced.string, collapse = ''), split = ' '))
  
  if(length(check.subs) > 1) {
    for(j in 2:(length(check.subs))) {
      
      split.entry <- unlist(strsplit(check.subs[j], split = ''))
      reduced.entry <- NULL # Store the column names here
      for(i2 in 1:(length(split.entry)-1)) {
        
        if(split.entry[i2] != " ") reduced.entry[length(reduced.entry)+1] <- split.entry[i2]
        if(i2 > 1 & split.entry[i2] != " " & split.entry[i2+1] == " ") reduced.entry[length(reduced.entry)+1] <- " "
        if(i2 == (length(split.entry)-1)) reduced.entry[length(reduced.entry)+1] <- split.entry[i2+1]
        
      }
      
      SLURM_OUT$JOBID[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[1]
      SLURM_OUT$PARTITION[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[2]
      SLURM_OUT$NAME[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[3]
      SLURM_OUT$USER[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[4]
      SLURM_OUT$ST[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[5]
      SLURM_OUT$TIME[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[6]
      SLURM_OUT$NODES[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[7]
      SLURM_OUT$`NODELIST(REASON)`[j-1] <- unlist(strsplit(paste(reduced.entry, collapse = ''), split = ' '))[8]
      
    }
  }
  
  return(SLURM_OUT)
  
}