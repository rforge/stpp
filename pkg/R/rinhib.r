rinhib <- function(npoints,s.region,t.region,hs="step",ps="min",thetas=0,deltas,ht="step",pt="min",thetat=1,deltat,recent="all",nsim=1,discrete.time=FALSE,replace=FALSE,inhibition=TRUE,...)
  {
  #
  # Simulate an inhibition or a contagious point process in a region D x T.
  #
  # Requires Splancs.
  #  
  # Arguments:
  #  
  #  s.region: two columns matrix specifying polygonal region containing
  #        all data locations. If poly is missing, the unit square is
  #        considered.
  #
  #  t.region: vector containing the minimum and maximum values of
  #            the time interval.
  #
  #         h: a function of the distance between points and theta.
  #            If inhibition=TRUE, h is monotone, increasing, and must tend
  #            to 1 when the distance tends to infinity. 0<= h(d,theta) <= 1.
  #            Else, h is monotone, decreasing, and must tend
  #            to 1 when the distance tends to 0. 0<= h(d,theta) <= 1.
  #  
  #         p: a function among "min", "max", "prod".  
  #  
  #   replace: logical allowing times repetition.
  #
  #   npoints: number of points to simulate. 
  #
  #      nsim: number of simulations to generate. Default is 1.
  #
  #
  # Value:
  #  pattern: list containing the points (x,y,t) of the simulated point
  #           process.
  #
  ##
  ## E. GABRIEL, 09/03/2007
  ##
  ## last modification: 13/03/2007
  ##                    27/03/2007   
  ##
  ##
  
  if (is.null(npoints)) stop("please specify the number of points to generate")
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  if (!(is.function(hs)))
    {
      if ((inhibition==TRUE) && (thetas==0) && (hs=="step") && (npoints * pi * deltas^2/4 > areapl(s.region)))
        stop(paste("s.region is too small to fit", npoints, "points", "at minimum distance", delta[1]))
    
      if ((inhibition==TRUE) && (thetas==0) && (hs=="step") && ((max(t.region)-min(t.region))/deltat<npoints))
        stop(paste("t.region is too small to fit", npoints, "points", "at minimum time interval", deltat))
    }
    
  pattern <- list()
  ni <- 1
  while(ni<=nsim)
    {
      pattern.interm <- spatial.inhibition(npoints,h=hs,theta=thetas,delta=deltas,p=ps,recent=recent,s.region=s.region,inhibition=inhibition,...)$pts
      times.interm <- temporal.inhibition(npoints,h=ht,theta=thetat,delta=deltat,p=pt,recent=recent,t.region=t.region,inhibition=inhibition,discrete.time=discrete.time,replace=replace,...)$times
      
      if (nsim==1)
        {
          pattern <- cbind(pattern.interm,times.interm)
          ni <-  ni+1
        }
      else
        {
          pattern[[ni]] <- cbind(pattern.interm,times.interm)
          ni <- ni+1
        }
    }
  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region)))
}




