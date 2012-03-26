gauss3D <- function(nx=100,ny=100,nt=100,xlim,ylim,tlim,separable=TRUE,model="exponential",param=c(1,1,1,1,1,2),scale=c(1,1),var.grf=1,mean.grf=0,exact=TRUE)
{

  N <- c(nx,ny,nt)

  mod <- set.cov(separable,model,param,var.grf)
  
  g <- floor(log(2*(N-1))/log(2))+1
  M <- 2^g

  count <- 1
  changG <- TRUE
  while(changG==TRUE)
    {
      L <- rep(-9999,M[1]*M[2]*M[3])  
      
      storage.mode(L) <- "double"
      
      res <- .Fortran("circ",
                      (L),
                      as.integer(M),
                      as.integer(N),
                      as.double(xlim),
                      as.double(ylim),
                      as.double(tlim),
                      as.integer(mod), 
                      as.double(param),
                      as.double(var.grf),
                      as.double(scale))[[1]]
 
      L <- array(res,dim=M)
      FTL <- fft(L)

      if (isTRUE(exact))
        {      
          if (min(Re(FTL))<0)
            {
              g <- g+1
              M <- 2^g
              changG <- TRUE
              count <- count+1
            }
          else
            changG <- FALSE
        }
      else
        {
          FTL[Re(FTL)<0] <- 0
          changG <- FALSE
        }

    }
  print(count)
  
  X <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
  X <- fft(X,inverse=TRUE)
  A <- sqrt(FTL)*X
  G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
  G <- G+mean.grf

## Remarks
# --------
#
#  The right way is:
#    X1 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#    X2 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#    X <- fft(X1+1i*X2,inverse=TRUE)
#    A <- sqrt(FTL)*X
#    G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
#  but it doesn't change the results and it is slightly faster. 
#  
# Taking the first elements of FTL and then computing G (changing M by N)
# is faster than taking the first elements of G, but it does not provide
# the same results
  
  return(G)
}



rlgcp <- function(s.region, t.region, replace=TRUE, npoints=NULL, nsim=1, nx=100, ny=100, nt=100,separable=TRUE,model="exponential",param=c(1,1,1,1,1,2),scale=c(1,1),var.grf=1,mean.grf=0,lmax=NULL,discrete.time=FALSE,exact=FALSE)
{
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)
  
  t.region <- sort(t.region)
  s.area <- areapl(s.region)
  t.area <- t.region[2]-t.region[1]
  tau <- c(start=t.region[1],end=t.region[2],step=(t.region[2]-t.region[1])/(nt-1))
  bpoly <- bbox(s.region)

  lambdamax <- lmax
  pattern <- list()
  index.t <- list()
  Lambdafin <- list()
  ni <- 1

  s.grid <- make.grid(nx,ny,s.region)
  s.grid$mask <- matrix(as.logical(s.grid$mask),nx,ny)

  if (discrete.time==TRUE)
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

  while(ni<=nsim)
    {
      S <- gauss3D(nx=nx,ny=ny,nt=nt,xlim=range(s.region[,1]),ylim=range(s.region[,2]),tlim=range(t.region),separable=separable,model=model,param=param,scale=scale,var.grf=var.grf,mean.grf=mean.grf,exact=exact)

      Lambda <- exp(S)

      mut <- rep(0,nt)
      for (it in 1:nt)
        {
          Lambda[,,it][s.grid$mask==FALSE] <- NaN
          mut[it] <- sum(Lambda[,,it],na.rm=TRUE)
        }
      
      if (is.null(npoints))
        {
          en <- sum(Lambda,na.rm=TRUE)*s.grid$xinc*s.grid$yinc
          npoints <- round(rpois(n=1,lambda=en),0)
        }

      if (is.null(lambdamax))
        lambdamax <- max(Lambda,na.rm=TRUE)
  
      npts <- round(lambdamax/(s.area*t.area),0)
      if (npts==0) stop("there is no data to thin")
  
      if ((replace==FALSE) & (nt < max(npts,npoints))) stop("when replace=FALSE, nt must be greater than the number of points used for thinning")

      if (discrete.time==TRUE)
        {
          vect <- seq(floor(t.region[1]),ceiling(t.region[2]),by=1)
          times.init <- sample(vect,nt,replace=replace)
        }
      else
        times.init <- runif(nt,min=t.region[1],max=t.region[2])
      
      samp <- sample(1:nt,npts,replace=replace,prob=mut/max(mut,na.rm=TRUE))
      times <- times.init[samp]

      retain.eq.F <- FALSE
      while(retain.eq.F==FALSE)
        {
          xy <- matrix(csr(poly=s.region,npoints=npts),ncol=2)
          x <- xy[,1]
          y <- xy[,2]
      
          prob <- NULL
          for(ix in 1:length(x))
            {
              nix <- iplace(X=s.grid$x,x=x[ix],xinc=s.grid$xinc)
              niy <- iplace(X=s.grid$y,x=y[ix],xinc=s.grid$yinc)
              nit <- iplace(X=t.grid$times,x=times[ix],xinc=t.grid$tinc)
              prob <- c(prob,Lambda[nix,niy,nit]/lambdamax)
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
          if (sum(retain==FALSE)==length(retain)) retain.eq.F <- FALSE
          else retain.eq.F <- TRUE
        }
      
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
          if (dim(xy)[2]==1) {wx <- xy[1]; wy <- xy[2]}
          else {wx <- xy[,1]; wy <- xy[,2]}
          
          if (isTRUE(replace))
            wsamp <- sample(1:nt,npoints-neffec,replace=replace,prob=mut/max(mut,na.rm=TRUE))
          else
            wsamp <- sample(samp.remain,npoints-neffec,replace=replace,prob=mut[samp.remain]/max(mut[samp.remain],na.rm=TRUE))
          
          wtimes <- times.init[wsamp]
          prob <- NULL
          for(ix in 1:length(wx))
            {
              nix <- iplace(X=s.grid$x,x=wx[ix],xinc=s.grid$xinc)
              niy <- iplace(X=s.grid$y,x=wy[ix],xinc=s.grid$yinc)
              nit <- iplace(X=t.grid$times,x=wtimes[ix],xinc=t.grid$tinc)
              prob <- c(prob,Lambda[nix,niy,nit]/lambdamax)
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
      pattern.interm <- cbind(x=x,y=y,t=times)

      if (nsim==1)
        {
          pattern <- as.3dpoints(pattern.interm)
          index.t <- index.times
          Lambdafin <- Lambda
        }
      else
        {
          pattern[[ni]] <- as.3dpoints(pattern.interm)
          index.t[[ni]] <- index.times
          Lambdafin[[ni]] <- Lambda
        }
      ni <- ni+1
    }

  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region,Lambda=Lambdafin,index.t=index.t)))
}









