iplace <- function(X,x,xinc)
  {
    n <- length(X)
    i <- 0
     repeat
       {
         i <- i+1
         if ((x >= X[i]-xinc) & (x < X[i]+xinc))
           break
       }
    ip <- i
    return(ip)
  }

