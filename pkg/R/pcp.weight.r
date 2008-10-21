pcp.weight <- function(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, mc=NULL, nsim=1, cluster="uniform", maxrad, infecD=TRUE, path=myfortranpath, ...)
{
  if (missing(cluster)) cluster <- "uniform"
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  if (missing(maxrad)) maxrad <- c(0.05,0.05)
  maxrads <- maxrad[1]
  maxradt <- maxrad[2]
  
  t.region <- sort(t.region)
  s.area <- areapl(s.region)
  t.area <- t.region[2]-t.region[1]
  
  if (is.null(nparents))
    {
      if (is.null(lambda)) stop("please specify either the number of parents or the intensity of parents process")
      if (is.function(lambda)) stop("please specify the number of parents")
      npar <- rpois(1,lambda)
    }
  else npar <- nparents
  if (is.null(lambda)) lambda <- npar/(s.area*t.area)
  
  parpts <- rpp(lambda=lambda,s.region=s.region,t.region=t.region,npoints=npar,...)$pts
  npoly <- length(s.region[, 1])
  polyx <- c(s.region[, 1], s.region[1, 1])
  polyy <- c(s.region[, 2], s.region[1, 2])
  binft <- t.region[1]
  bsupt <- t.region[2]
  W <- rep(0,npar)
  cont <- as.numeric(infecD)
#  dyn.load("/home/gabriel/functions/SimulSTPP/libF/edgecorrect.dll")
 pathname <- paste(path,"edgecorrect.dll",sep="")
  dyn.load(pathname)
  weight <- .Fortran("edgecorrect", as.double(parpts[,1]),
                     as.double(parpts[,2]), as.double(parpts[,3]), 
                     as.integer(npar), as.double(polyx),
                     as.double(polyy), as.integer(npoly),
                     as.double(maxrads), as.double(maxradt),
                     as.double(binft), as.double(bsupt),
                     as.integer(cont), as.double(W))
  W <- weight[[13]]

  if (is.null(npoints))
    {
      if (is.null(mc)) stop("please specify either the number of points to simulate or the mean number of children per parents")
      npts <- round(sum(rpois(npar,lambda=mc)*W),0)
    }
  else npts <- npoints
  pattern.interm <- matrix(0,npts,3)
  nchild <- matrix(0,npar,1)
  ipt <- 1
  while(ipt <= npts)
    {
      ipar <- sample(1:npar,1,prob=W)
      xpar <- parpts[ipar,1]
      ypar <- parpts[ipar,2]
      zpar <- parpts[ipar,3]
      nchild[ipar,1] <- nchild[ipar,1]+1
      
      if (length(cluster)==1)
        {
          s.distr <- cluster
          t.distr <- cluster
        }
      else
        {
          s.distr <- cluster[1]
          t.distr <- cluster[2]
        }
      outcirc <- 1
      while(outcirc == 1)
        {
          if (s.distr=="uniform")
            {
              xp <- xpar+runif(1,min=-maxrads,max=maxrads)
              yp <- ypar+runif(1,min=-maxrads,max=maxrads)
            }
          if (s.distr=="normal")
            {
                  xp <- xpar+rnorm(1,mean=0,sd=maxrads/2)
                  yp <- ypar+rnorm(1,mean=0,sd=maxrads/2)
                }
          if (s.distr=="exponential")
            {
              xp <- xpar+rexp(1,rate=1/maxrads)
              yp <- ypar+rexp(1,rate=1/maxrads)
            }
          if (t.distr=="uniform")
            {
              if (infecD==TRUE) 
                zp <- zpar+runif(1,min=-maxradt,max=maxradt)
              else
                zp <- zpar+runif(1,min=0,max=maxradt)
            }
          if (t.distr=="normal")
            {
              if (infecD==TRUE) 
                zp <- zpar+abs(rnorm(1,mean=0,sd=maxradt/2))
              else
                zp <- zpar+rnorm(1,mean=0,sd=maxradt/2)
            }
          if (t.distr=="exponential")
            {
              if (infecD==TRUE) 
                zp <- zpar+rexp(1,rate=1/maxradt)
              else
                zp <- zpar+sample(c(-1,1),1)*rexp(1,rate=1/maxradt)
            }
              
          if ((inout(as.points(x=xp,y=yp),poly=s.region)==TRUE) & (zp > t.region[1] & zp < t.region[2]) & (sqrt( (((xpar-xp)^2)/maxrads^2) + (((ypar-yp)^2)/maxrads^2)) < 1)) outcirc <- 0
        }
      pattern.interm[ipt,1:3] <- cbind(xp,yp,zp)
      ipt <- ipt+1
    }
  pts <- pattern.interm
}
  
