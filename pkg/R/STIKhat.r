STIKhat <- function(xyt, s.region, t.region, dist, times, lambda, correction = TRUE, infectious = TRUE) 
{

#  library(splancs)
  #
  # This function computes the space-time inhomogeneous K function.
  #
  # Arguments:
  #  - xyt: coordinates and times of the point process,
  #  - s.region: set of points defining the domain D,
  #  - t.region: vector giving the lower and upper boundaries of T
  #  - dist: vector of distances h at which K() is computed,
  #  - times: vector of times t at which K() is computed,
  #  - lambda: vector of values of the space-time intensity
  #            function evaluated at the points of the pattern,
  #  - correction: logical value. If True, the Ripley's edge correction
  #                is used.
  #  - infectious: logical value. If True, the STIK function is the one defined
  #            in the case of infectious diseases.
  #
  # Value:
  #  - Khat: a ndist x ntimes matrix containing the values of K(h,t),
  #  - dist: the vector of distances used,
  #  - times: the vector of times used.
  #
  #
  # E. Gabriel, september 2005
  #
  # last modification: 25.10.2006
  #

  if (missing(s.region)) s.region <- sbox(xyt[,1:2],xfrac=0.01,yfrac=0.01)
  if (missing(t.region)) t.region <- range(xyt[,3],na.rm=T)
    
  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[, 1]
  ptsy <- pts[, 2]
  ptst <- xytimes
  npt <- length(ptsx)
  ndist <- length(dist)
  dist <- sort(dist)
  ntimes <- length(times)
  times <- sort(times)
  bsupt <- max(t.region)
  binft <- min(t.region)
    
  area <- areapl(s.region)*(bsupt-binft)

  np <- length(s.region[, 1])
  polyx <- c(s.region[, 1], s.region[1, 1])
  polyy <- c(s.region[, 2], s.region[1, 2])
  hkhat <- array(0, dim = c(ndist,ntimes))

  if(missing(lambda))
    {
      misl <- 1
      lambda <- rep(1,npt)
    }
  else misl <- 0
  if (length(lambda)==1) lambda <- rep(lambda,npt)
  if (correction==T){edg <- 1} else  edg <- 0
  if (infectious==T){infd <- 1} else infd <- 0
  storage.mode(hkhat) <- "double"

  if (area>10) lambda <- lambda*area
  
#  dyn.load("/home/gabriel/functions/STIKfunction/libF/stikfunction.dll")
#  pathname <- paste(path,"stikfunction.dll",sep="")
#  dyn.load(pathname)
  nev <- rep(0,ntimes)
  klist <- .Fortran("stikfunction", as.double(ptsx),
                    as.double(ptsy), as.double(ptst), 
                    as.integer(npt), as.double(polyx),
                    as.double(polyy), as.integer(np),
                    as.double(dist), as.integer(ndist),
                    as.double(times), as.integer(ntimes),
                    as.double(bsupt), as.double(binft),
                    as.double(lambda), as.integer(infd),
                    as.integer(edg), (hkhat))
  khat <- klist[[17]]
  if (misl==1) 
   {
    if (area>10)
        hkhat <- ((area^3)/(npt*(npt-1)))*khat
      else
        hkhat <- (area/(npt*(npt-1)))*khat
    }
  else
    {
      hkhat <- khat
        
      if (area>10)
        hkhat <- hkhat*area
      else
        hkhat <- hkhat/area
    }
    if(dist[1]==0) hkhat[1,]=0
    if(times[1]==0) hkhat[,1]=0
  
  Khpp <- matrix(0,ncol=length(times),nrow=length(dist))
  for(i in 1:length(dist)){Khpp[i,]<-pi*(dist[i]^2)*times}
  
  if (infectious==T) Ktheo <- Khpp
  else Ktheo <- 2*Khpp
 
  lhat <- hkhat-Ktheo

  invisible(return(list(Khat=hkhat,Ktheo=Ktheo,dist=dist,times=times,infectious=infectious)))  
}
