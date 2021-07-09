in.polygon <- function(points = NULL, polygon = NULL) {

  names(polygon) <- names(points) <- c('x', 'y')
  
  closure.test <-
    polygon[1,1] == polygon[nrow(polygon),1] &
    polygon[1,2] == polygon[nrow(polygon),2]
  
  if(!closure.test & nrow(polygon) > 2) {
    message('Passed object not a polygon or is not closed.')
    return(NA)
  }
  
  for(i in 1:nrow(points)) {
    
    test.point <- points[i,]
      
      count <- 0
      for(j in 1:(nrow(polygon)-1)) {
        
        p.point.1 <- polygon[j,]
        p.point.2 <- polygon[j+1,]
        
        p.min <- which.min(c(p.point.1$x, p.point.2$x))
        
        if(p.min == 1 & !(p.point.1$x == p.point.2$x)) {
          
          m = (p.point.2$y - p.point.1$y)/(p.point.2$x - p.point.1$x)
          b = p.point.1$y - m*p.point.1$x
          
          cross.x = (test.point$y - b)/m
          cross.y = m*cross.x + b
          
          if(m == 0 & (test.point$y == p.point.1$y) &
             test.point$x <= max(p.point.1$x, p.point.2$x) &
             test.point$x >= min(p.point.1$x, p.point.2$x)) {
            count <- count + 1
          } else if(m != 0) {
            if(cross.x >= test.point$x &
             cross.y <= max(p.point.1$y, p.point.2$y) &
             cross.y >= min(p.point.1$y, p.point.2$y)) {count <- count + 1}
          }
          
        } else if (p.min == 2 & !(p.point.1$x == p.point.2$x)) {
          
          m = (p.point.1$y - p.point.2$y)/(p.point.1$x - p.point.2$x)
          b = p.point.2$y - m*p.point.2$x
          
          cross.x = (test.point$y - b)/m
          cross.y = m*cross.x + b
          
          if(m == 0 & (test.point$y == p.point.1$y) &
             test.point$x <= max(p.point.1$x, p.point.2$x) &
             test.point$x >= min(p.point.1$x, p.point.2$x)) {
            count <- count + 1
          } else if(m != 0) {
            if(cross.x >= test.point$x &
             cross.y <= max(p.point.1$y, p.point.2$y) &
             cross.y >= min(p.point.1$y, p.point.2$y)) {count <- count + 1}
          }
          
        } else if (p.point.1$x == p.point.2$x) {
          
          if(test.point$y <= max(p.point.1$y, p.point.2$y) &
             test.point$y >= min(p.point.1$y, p.point.2$y) &
             test.point$x <= p.point.1$x) {count <- count + 1}
            
        }

      }
      
      # remove test points crossing at vertices
      crossed.vert <- subset(polygon,
                             y == test.point$y & x >= test.point$x)
      
      count <- count - nrow(crossed.vert)
      
      #check the count here. Odd or even?
      if((count %% 2) == 0 | count == 0) points$polygon[i] <- 'out'
      if((count %% 2) != 0) points$polygon[i] <- 'in'
      
  }
  
  return(points)
  
}
