#separate camelback notation
separate.camelback <- function(strings = NULL) {
  
  output <- NULL
  for(i in 1:length(strings)) {
    
    string <- strings[i]
    
    Upper.Idx <- unlist(gregexpr("[A-Z]", string))
    Upper.Idx <- Upper.Idx[Upper.Idx != 1]
    
    for(j in 1:length(Upper.Idx)) {
      
      if(j != 1) {
        Upper.Idx <- unlist(gregexpr("[A-Z]", string))
        Upper.Idx <- Upper.Idx[Upper.Idx != 1]
      }
      
      current <- substr(string, Upper.Idx[j]-1, Upper.Idx[j])
      replace <- paste(substr(current, 1, 1), substr(current, 2, 2))
      
      string <- gsub(current, replace, string)
      
    }
    
    output[i] <- string
    
  } #close string loop
  
  output
  
}