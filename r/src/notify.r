notify <- function(user.email = NULL, subject.text = NULL, body.text = NULL) {
  
  system(paste0("echo", " '", body.text, "' ", "|", " mail -s", " '", subject.text, "' ", user.email))
  
}



