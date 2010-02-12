rpp <- function(lambda, s.region, t.region, npoints=NULL, nsim=1, replace=TRUE, discrete.time=FALSE, nx=100, ny=100, nt=100, lmax=NULL, Lambda=NULL, mut=NULL, ...)
{
  #
  # Simulate a space-time Poisson process in a region D x T.
  #
  # Requires Splancs package.
  #  
  # Arguments:
  #
  #        lambda: Spatio-temporal intensity of the Poisson process. 
  #                If 'lambda' is a single positive number, the function 
  #			 generates realisations of a homogeneous Poisson process, 
  # 			 whilst if 'lambda' is a function of the form 
  #			 lambda(x,y,t,...) or a character it generates 
  #   		 realisations of an inhomogeneous Poisson process.
  #
  #      s.region: two columns matrix specifying polygonal region containing
  #                all data locations. If s.region is missing, the unit square
  #                is considered.
  #
  #      t.region: vector containing the minimum and maximum values of
  #                the time interval. If t.region is missing, the interval
  #                [0,1] is considered.
  #
  #       replace: logical allowing times repetition.
  #
  #       npoints: number of points to simulate. If NULL (default), the
  #                number of points is from a Poisson distribution with
  #                mean the double integral of lambda over s.region and
  #                t.region.
  #
  # discrete.time: if TRUE, times belong to N, otherwise belong to R+.
  #
  #          nsim: number of simulations to generate. Default is 1.
  #
  #      nx,ny,nt: define the 3-D grid on which the intensity is evaluated.
  #
  #          lmax: upper bound for the value of lambda(x,y,t), if lambda
  #                is a function.
  #
  #	      Lambda: matrix of spatial intensity if 'lambda' is a character.
  #
  # 		   mut: vector of temporal intensity if 'lambda' is a character.
  #
  #
  # Value:
  #     xyt: matrix (or list if nsim>1) containing the points (x,y,t)
  #           of the simulated point process.
  #  Lambda: an array of the intensity surface at each time.
  #
  # NB: the probability that an event occurs at a location s and a time t
  #     is:
  #     p(s,t)=lambda(s,t)/(max_{(s_i,t_i), i=1,...,ngrid} lambda(s_i,t_i)).
  #

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

  lambdamax <- lmax
  pattern <- list()
  index.t <- list()
  ni <- 1

  #
  # Homogeneous Poisson Process
  #

  if (is.numeric(lambda))
    {
      while(ni<=nsim)
        {
          hpp <- rhpp(lambda=lambda, s.region=s.region, t.region=t.region, npoints=npoints, replace=replace, discrete.time=discrete.time)
          if (nsim==1)
            {
              pattern <- as.3dpoints(hpp$pts)
              index.t <- hpp$index.t
            }
          else
            {
              pattern[[ni]] <- as.3dpoints(hpp$pts)
              index.t[[ni]] <- hpp$index.t
            }
          ni <- ni+1
        }
      Lambda <- NULL
    }
    
  #
  # Inhomogeneous Poisson Process
  #

  if (is.function(lambda))
    {
      s.grid <- make.grid(nx,ny,s.region)
      s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)
      if (discrete.time==T)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          if (nt>length(vect))
            {
              nt <- length(vect)
              warning("nt used is less than the one given in argument")
              t.grid <- list(times=vect,tinc=1)
            }
          else
            {
              vect <- round(seq(floor(t.region[1]),ceiling(t.region[2]),length=nt))
              t.grid <- list(times=vect,tinc=round(t.area/(nt-1)))
            }
        }
      else
        t.grid <- list(times=seq(t.region[1],t.region[2],length=nt),tinc=(t.area/(nt-1)))

      Lambda <- array(NaN,dim=c(nx,ny,nt)) 
      mut <- rep(0,nt)
 #         maxlambda <- rep(0,nt)
      for(it in 1:nt)
        {
          L <- lambda(as.vector(s.grid$X),as.vector(s.grid$Y),t.grid$times[it],...)
          M <- matrix(L,ncol=ny,nrow=nx,byrow=T)
          M[!(s.grid$mask)] <- NaN
          Lambda[,,it] <- M
 #             maxlambda[it] <- max(Lambda[,,it],na.rm=T)
          mut[it] <- sum(Lambda[,,it],na.rm=T)
        }
          
      while(ni<=nsim)
        {
          ipp <- ripp(lambda=lambda, s.region=s.region, t.region=t.region, npoints=npoints, replace=replace, discrete.time=discrete.time, nx=nx, ny=ny, nt=nt, lmax=lmax, Lambda=Lambda, mut=mut, ...)
          
          if (nsim==1)
            {
              pattern <- as.3dpoints(ipp$pts)
              index.t <- ipp$index.t
            }
          else
            {
              pattern[[ni]] <- as.3dpoints(ipp$pts)
              index.t[[ni]] <- ipp$index.t
            }
          ni <- ni+1
        }
    }

  if (is.character(lambda))
    {
      while(ni<=nsim)
        {
          ipp <- ripp(lambda=lambda, s.region=s.region, t.region=t.region, npoints=npoints, replace=replace, discrete.time=discrete.time, nx=nx, ny=ny, nt=nt, lmax=lmax, Lambda=Lambda, mut=mut, ...)
          
          if (nsim==1)
            {
              pattern <- as.3dpoints(ipp$pts)
              index.t <- ipp$index.t
            }
          else
            {
              pattern[[ni]] <- as.3dpoints(ipp$pts)
              index.t[[ni]] <- ipp$index.t
            }
          ni <- ni+1
        }
    }
  
  invisible(return(list(xyt=pattern,index.t=index.t,s.region=s.region,t.region=t.region,lambda=lambda,Lambda=Lambda)))
}



