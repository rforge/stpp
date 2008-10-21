ripp <- function(lambda, s.region, t.region, npoints=NULL, replace=T, discrete.time=FALSE, nx=100, ny=100, nt=100, lmax=NULL, Lambda=NULL, mut=NULL, ...)
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

  lambdamax <- lmax
  pattern <- list()
  index.t <- list()

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

      if (is.null(Lambda) & is.null(mut))
        {
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
        }
      
      if (is.null(npoints)==T)
        {
          if (t.area==0)
            { en <- sum(Lambda,na.rm=T)*s.grid$xinc*s.grid$yinc }
          else
            {
              en <- sum(Lambda,na.rm=T)*s.grid$xinc*s.grid$yinc*t.grid$tinc 
              npoints <- round(rpois(n=1,lambda=en),0)
            }
        }

      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=T)
      npts <- round(lambdamax/(s.area*t.area),0)
      if (npts==0) stop("there is no data to thin")
      
      xy <- matrix(csr(poly=s.region,npoints=npts),ncol=2)
      x <- xy[,1]
      y <- xy[,2]
      
      if ((replace==F) && (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")
      if (discrete.time==T)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
        }
      else
        times.init <- runif(nt,min=t.region[1],max=t.region[2])
      
      samp <- sample(1:nt,npts,replace=replace,prob=mut/max(mut,na.rm=T))
      times <- times.init[samp]
      prob <-  lambda(x,y,times,...)/lambdamax
      u <- runif(npts)
      retain <- u <= prob
      if (sum(retain==F)==length(retain))
        {
          lambdas <- matrix(0,nrow=nx,ncol=ny)
          for(ix in 1:nx){for(iy in 1:ny){
            lambdas[ix,iy] <- median(Lambda[ix,iy,],na.rm=T)}}
          lambdamax <- max(lambdas,na.rm=T)
          prob <-  lambda(x,y,times,...)/lambdamax
          retain <- u <= prob
          if (sum(retain==F)==length(retain)) stop ("no point was retained at the first iteration, please check your parameters")
        }
      x <- x[retain]
      y <- y[retain]
      samp <- samp[retain]
      samp.remain <- (1:nt)[-samp]
      times <- times[retain]
      
      neffec <- length(x)
      while(neffec < npoints)
        {
          xy <- as.matrix(csr(poly=s.region,npoints=npoints-neffec))
          if(dim(xy)[2]==1){wx <- xy[1]; wy <- xy[2]}
          else{wx <- xy[,1]; wy <- xy[,2]}
          if(replace==F)
            { wsamp <- sample(samp.remain,npoints-neffec,replace=replace,prob=mut[samp.remain]/max(mut[samp.remain],na.rm=T)) }
          else{ wsamp <- sample(1:nt,npoints-neffec,replace=replace,prob=mut/max(mut,na.rm=T)) }
          wtimes <- times.init[wsamp]
#              lambdamax <- maxlambda[wsamp]
          prob <-  lambda(wx,wy,wtimes,...)/lambdamax
          u <- runif(npoints-neffec)
          retain <- u <= prob
          x <- c(x,wx[retain])
          y <- c(y,wy[retain])
          times <- c(times,wtimes[retain])
          samp <- c(samp,wsamp[retain])
          samp.remain <- (1:nt)[-samp]
          neffec <- length(x)
        }
      times <- sort(times)
      index.times <- sort(samp)
      pattern.interm <- list(x=x,y=y,t=times,index.t=index.times)
          pattern <- cbind(x=pattern.interm$x,y=pattern.interm$y,t=pattern.interm$t)
      index.t <- pattern.interm$index.t
    }

  if (is.character(lambda))
    {
      if (is.null(Lambda) & is.null(mut))
        stop("Lambda and mut must be specified")
      nx <- dim(Lambda)[1]
      ny <- dim(Lambda)[2]
      nt <- length(mut)
      
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
          
      if (is.null(npoints))
        {
          en <- sum(Lambda,na.rm=T)*s.grid$xinc*s.grid$yinc
          npoints <- round(rpois(n=1,lambda=en),0)
        }
      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=T)*max(mut)/npoints
#      npts <- round(lambdamax/(s.area*t.area),0)
      npts <- npoints
      if (npts==0) stop("there is no data to thin")
     
      if ((replace==F) && (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")
      if (discrete.time==T)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
        }
      else
        times.init <- runif(nt,min=t.region[1],max=t.region[2])
      
      samp <- sample(1:nt,npts,replace=replace,prob=mut/max(mut,na.rm=T))
      times <- times.init[samp]

      retain.eq.F <- F
      while(retain.eq.F==F)
        {
          xy <- matrix(csr(poly=s.region,npoints=npts),ncol=2)
          x <- xy[,1]
          y <- xy[,2]
          
          prob <- NULL
          for(nx in 1:length(x))
            {
              nix <- iplace(X=s.grid$x,x=x[nx],xinc=s.grid$xinc)
              niy <- iplace(X=s.grid$y,x=y[nx],xinc=s.grid$yinc)
              nit <- iplace(X=t.grid$times,x=times[nx],xinc=t.grid$tinc)
              prob <- c(prob,(Lambda[nix,niy]*mut[nit]/npoints)/lambdamax)
            }
          
          M <- which(is.na(prob))
          if (length(M)!=0)
            {
              x <- x[-M]
              y <- y[-M]
              times <- times[-M]
              prob <- prob[-M]
              npts <- length(x)
            }
          
          u <- runif(npts)
          retain <- u <= prob
          if (sum(retain==F)==length(retain)) retain.eq.F <- F
          else retain.eq.F <- T
        }
#      if (sum(retain==F)==length(retain)) stop ("no point was retained at the first iteration, please check your parameters")
      
      x <- x[retain]
      y <- y[retain]
      samp <- samp[retain]
      samp.remain <- (1:nt)[-samp]
      times <- times[retain]

      neffec <- length(x)
      if (neffec > npoints)
        {
          retain <- 1:npoints
          x <- x[retain]
          y <- y[retain]
          samp <- samp[retain]
          samp.remain <- (1:nt)[-samp]
          times <- times[retain]
        }

      while(neffec < npoints)
        {
          xy <- as.matrix(csr(poly=s.region,npoints=npoints-neffec))
          if(dim(xy)[2]==1){wx <- xy[1]; wy <- xy[2]}
          else{wx <- xy[,1]; wy <- xy[,2]}
          if(replace==F)
            { wsamp <- sample(samp.remain,npoints-neffec,replace=replace,prob=mut[samp.remain]/max(mut[samp.remain],na.rm=T)) }
          else{ wsamp <- sample(1:nt,npoints-neffec,replace=replace,prob=mut/max(mut,na.rm=T)) }
          wtimes <- times.init[wsamp]
#              lambdamax <- maxlambda[wsamp]

          prob <- NULL
          for(nx in 1:length(wx))
            {
              nix <- iplace(X=s.grid$x,x=wx[nx],xinc=s.grid$xinc)
              niy <- iplace(X=s.grid$y,x=wy[nx],xinc=s.grid$yinc)
              nit <- iplace(X=t.grid$times,x=wtimes[nx],xinc=t.grid$tinc)
              prob <- c(prob,(Lambda[nix,niy]*mut[nit]/npoints)/lambdamax)
            }
          M <- which(is.na(prob))
          if (length(M)!=0)
            {
              wx <- wx[-M]
              wy <- wy[-M]
              wtimes <- wtimes[-M]
              prob <- prob[-M]
            }
          if (neffec > 0)
            {
              u <- runif(length(prob))
              retain <- u <= prob
              x <- c(x,wx[retain])
              y <- c(y,wy[retain])
              times <- c(times,wtimes[retain])
              samp <- c(samp,wsamp[retain])
              samp.remain <- (1:nt)[-samp]
              neffec <- length(x)
            }
        }
      times <- sort(times)
      index.times <- sort(samp)
      pattern.interm <- list(x=x,y=y,t=times,index.t=index.times)
      pattern <- cbind(x=pattern.interm$x,y=pattern.interm$y,t=pattern.interm$t)
      index.t <- pattern.interm$index.t
    }
      
  invisible(return(list(pts=pattern,index.t=index.t)))
}
  
