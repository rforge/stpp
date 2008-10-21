rhpp <- function(lambda, s.region, t.region, npoints=NULL, replace=TRUE, discrete.time=FALSE)
  {
    if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
    if (missing(t.region)) t.region <- c(0,1)
    
    s.area <- areapl(s.region)
    t.region <- sort(t.region)
    t.area <- t.region[2]-t.region[1]

    if (missing(lambda) & !(is.null(npoints)))
      {
        if (t.area==0) lambda <- npoints/(s.area)
        else lambda <- npoints/(s.area * t.area)
      }
    
    pattern <- list()
    index.t <- list()
   if (is.numeric(lambda))
    {
      if (is.null(npoints)==TRUE)
        {
          if (t.area==0)
            { npoints <- round(rpois(n=1,lambda=lambda * s.area),0) }
          else
            { npoints <- round(rpois(n=1,lambda=lambda * s.area * t.area),0) }
        }
      xy <- matrix(csr(poly=s.region,npoints=npoints),ncol=2)
      x <- xy[,1]
      y <- xy[,2]
      npoints <- length(x)
      
      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          if ((length(vect)<npoints) & (replace==F))
            stop("when replace=FALSE and discrete.time=TRUE, the length of seq(t.region[1],t.region[2],by=1) must be greater than the number of points")
          names(vect) <- 1:length(vect) 
          M <- sample(vect,npoints,replace=replace)
          times <- M
          names(times) <- NULL
          samp <- as.numeric(names(M))
        }
      else
        {
          times <- runif(npoints,min=t.region[1],max=t.region[2])
          samp <- sample(1:npoints,npoints,replace=replace)
          times <- times[samp]
        }
      times <- sort(times)
      index.times <- sort(samp)
      pattern.interm <- list(x=x,y=y,t=times,index.t=index.times)      
      pattern <- cbind(x=pattern.interm$x,y=pattern.interm$y,t=pattern.interm$t)
      index.t <- pattern.interm$index.t
    }
   else
     stop("lambda must be numeric")

    invisible(return(list(pts=pattern,index.t=index.t)))
  }
