SLURM.wait <- function(selected.partition = NULL, user.id = NULL) {
  #' monitor parallel jobs. A new job will not be submitted until the
  #' previous one is completed.
  user.jobs <- subset(SLURM.jobs(partition = selected.partition),
                      USER == user.id)
  msg.flag <- 0 # msg.flag == 0 will display a message when jobs are running.
  
  # Hang out while previous jobs are running.
  while(nrow(user.jobs) > 1) {
    if(msg.flag == 0) message('User has jobs in progress. Please wait.')
    Sys.sleep(300); msg.flag <- 1
    user.jobs <- subset(SLURM.jobs(partition = selected.partition),
                        USER == user.id)
  }; msg.flag <- 0
}