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
  #            to 1 when the distance tends to infinity. 0<= h(d,theta,delta) <= 1.
  #            Else, h is monotone, decreasing, and must tend
  #            to 1 when the distance tends to 0. 0<= h(d,theta,delta) <= 1.
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

  if (!(is.function(h)))
    {
      models <- c("step","gaussian")
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
		  if (inhibition==TRUE) res[d<=delta] <- theta
		  else res[d>=delta] <- theta
              return(res)
            }
        }
      if (h=="gaussian")
        {
          hk <- function(d,theta,delta)
            {
              if (inhibition==TRUE) 
			{
			res=NULL
			for(i in 1:length(d))
				{	
				if (d[i]<=delta) res=c(res,0)
				if (d[i]>(delta+theta/2)) res=c(res,1)
				if (d[i]>delta & d[i]<=(delta+theta/2)) res=c(res,exp(-((d[i]-delta-theta/2)^2)/(2*(theta/8)^2)))
				}
			}
		  else
			{
			res=NULL
			for(i in 1:length(d))
				{	
				if (d[i]<delta) res=c(res,1)
				else res=c(res,exp(-((d[i]-delta)^2)/(2*(theta/8)^2)))
				}
			}
	   	  return(res)
		}
	  }
	}
 else
   	{   
          hk <- function(d,theta,delta)
            {
            res <- h(d,theta,delta)
            return(res)
		}
	}

	pk <- function(d,h,recent,theta,delta)
        {
          if (recent=="all")
		{
		 if (p=="min") res <- min(h(d=d,theta=theta,delta=delta))
		 if (p=="max") res <- max(h(d=d,theta=theta,delta=delta))
		 if (p=="prod") res <- prod(h(d=d,theta=theta,delta=delta))
		}
          else
            {
              if (is.numeric(recent))
                {
                  if(recent<=length(d))
				{
				if (p=="min") res <- min(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				if (p=="max") res <- max(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				if (p=="prod") res <- prod(h(d=d[(length(d)-recent+1):length(d)],theta=theta,delta=delta))
				}
                  else
				{
				if (p=="min") res <- min(h(d=d,theta=theta,delta=delta))
				if (p=="max") res <- max(h(d=d,theta=theta,delta=delta))
				if (p=="prod") res <- prod(h(d=d,theta=theta,delta=delta))
				}
                }
              else
                stop("'recent' must be numeric")
            }
          return(res)
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
             umax <- pk(d=sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2),hk,recent,theta,delta)             
		
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
              if (all(sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2) < delta))
                umax <- 1            
              else		
                umax <- pk(d=sqrt((xy[1] - pattern.interm[,1])^2 + (xy[2] - pattern.interm[,2])^2),hk,recent,theta,delta)             	
                
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


