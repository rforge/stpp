spatial.inhibition <- function(npoints,h,theta,delta,p,recent="all",s.region,inhibition=TRUE,...)
  {
  #
  # Simulate an inhibition or a contagious spatial point process in a region D
  #
  # Requires Splancs.
  #  
  # Arguments:
  #  
  #  s.region: two columns matrix specifying polygonal region containing
  #        all data locations. If poly is missing, the unit square is
  #        considered.
  #
  #         h: a function of the distance between points.
  #            If inhibition=TRUE, h is monotone, increasing, and must tend
  #            to 1 when the distance tends to infinity. 0<= h(d,theta) <= 1.
  #            Else, h is monotone, decreasing, and must tend
  #            to 1 when the distance tends to 0. 0<= h(d,theta) <= 1.
  #  
  #         p: a function among "min", "max", "prod".  
  #  
  #   npoints: number of points to simulate. 
  #
  #
  # Value:
  #  pattern: list containing the points (x,y) of the simulated point
  #           process.
  #
  ##
  ## E. GABRIEL, 26/03/2007
  ##
  ## last modification: 21/10/2008
  ##
  ##
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)

  if ((inhibition==TRUE) && (npoints * pi * delta[1]^2/4 > areapl(s.region)))
        stop(paste("s.region is too small to fit", npoints, "points", "at minimum distance", delta[1]))
        
  if (!(is.function(h)))
    {
      models <- c("step","exponential","gaussian")
      if (sum(h==models)==0)
        {
          message <- paste("this model is not implemented, please choose among: ",paste(models,"",sep=" ",collapse="and "))
          stop(message)
        }
      if (h=="step")
        {
          hk <- function(d,theta,delta)
            {
              res <- rep(1,length(d))
              res[d<=delta] <- theta
              return(res)
            }
        }
      
      if (h=="exponential")
        {
          hk <- function(d,theta,delta)
            {
              res <- exp(d*theta)/max(exp(d*theta))
              return(res)
            }
        }
      if (h=="gaussian")
        {
          hk <- function(d,theta,delta)
            {
              res <- exp((d^2)*theta)/max(exp((d^2)*theta))
              return(res)
            }
        }
    }
 
  if (inhibition==FALSE)
    {
      if (!(is.function(h))) 
        {
          hh <- hk
          hk <- function(d,theta,delta)
            {
            res <- hh(max(d)-d,theta,max(d)-delta)
            return(res)
            }
        }  
      else 
        {
          hh <- h
          hk <- function(d,theta,delta)
            {
            res <- h(max(d)-d,theta,max(d)-delta,...)
            return(res)
            }
        }
    }

  if (p=="min")
    {
      pk <- function(d,h,recent,theta,delta,...)
        {
          if (recent=="all")
            res <- min(h(d=d,theta=theta,delta=delta,...))
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(d))
                    res <- min(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta,...))
                  else
                    res <- min(h(d=d,theta=theta,delta=delta,...))
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
        }
    }

    if (p=="max")
    {
      pk <- function(d,h,recent,theta,delta,...)
        {
          if (recent=="all")
            res <- max(h(d=d,theta=theta,delta=delta,...))
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(d))
                    res <- max(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta,...))
                  else
                    res <- max(h(d=d,theta=theta,delta=delta,...))
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
        }
    }

    if (p=="prod")
    {
      pk <- function(d,h,recent,theta,delta,...)
        {
          if (is.numeric(recent) && recent<=length(d))
            res <- prod(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta,...))
          else
            res <- prod(h(d=d,theta=theta,delta=delta,...))
          return(res)
        }
    }

  xy <- csr(n=1,poly=s.region)
  npts <- 1
  pattern.interm <- cbind(x=xy[1],y=xy[2])
  if (inhibition==TRUE)
    {
      while(npts < npoints)
        {
          prob <- runif(1)

          xy <- csr(n=1,poly=s.region)
          if (all((sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2)) > delta))
            umax <- 1
          else
            {
              if (is.function(h))
                umax <- pk(d=sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2),h,recent,theta,delta,...)
              else
                {
                  h <- hk
                  umax <- pk(d=sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2),h,recent,theta,delta,...)
                }
            }
          if (prob<umax)
            {
              pattern.interm <- rbind(pattern.interm,c(xy[1],xy[2]))
              npts <- npts+1
            }
        }
    }
  else
    {
      while(npts < npoints)
        {
          prob <- runif(1)

          continue <- FALSE
          while(continue==FALSE)
            {
              xy <- csr(n=1,poly=s.region)
              if (sqrt((xy[1] - pattern.interm[npts,1])^2 + (xy[2] - pattern.interm[npts,2])^2) < delta)
                umax <- 1            
              else
                {
                  if (is.function(h))
                    umax <- pk(d=sqrt((xy[1] - pattern.interm[npts,1])^2 + (xy[2] - pattern.interm[npts,2])^2),h,recent,theta,delta,...)
                  else
                    {
                      h <- hk
                      umax <- pk(d=sqrt((xy[1] - pattern.interm[npts,1])^2 + (xy[2] - pattern.interm[npts,2])^2),h,recent,theta,delta,...)
                    }
                }
              if (prob < umax)
                {
                  pattern.interm <- rbind(pattern.interm,c(xy[1],xy[2]))
                  npts <- npts+1
                  continue <- TRUE
                }
            }
        }
    }

  invisible(return(list(pts=pattern.interm,s.region=s.region)))
}


