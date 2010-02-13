pcp.larger.region <- function(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, mc=NULL, nsim=1, cluster="uniform", maxrad, infecD=TRUE, maxradlarger, ...)
{
  if (missing(cluster)) cluster <- "uniform"
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  if (missing(maxrad)) maxrad <- c(0.05,0.05)
  if (length(maxrad)==1) maxrad=rep(maxrad,2)
  maxrads <- maxrad[1]
  maxradt <- maxrad[2]

  if (missing(maxradlarger)) maxradlarger <- maxrad
  maxradlargers <- maxradlarger[1]
  maxradlargert <- maxradlarger[2]
  
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
  
  t.region <- sort(t.region)
  s.area <- areapl(s.region)
  t.area <- t.region[2]-t.region[1]

  if (maxradlargers==0)
    s.larger <- s.region
  else
    {
      s.larger <- rbind(s.region+maxradlargers,s.region-maxradlargers,cbind(s.region[,1]+maxradlargers,s.region[,2]-maxradlargers),cbind(s.region[,1]-maxradlargers,s.region[,2]+maxradlargers))
      M <- chull(s.larger)
      s.larger <- cbind(s.larger[M,1],s.larger[M,2])
    }
  t.larger <- c(t.region[1]-maxradlargert,t.region[2]+maxradlargert)

  if (t.larger[1]<0) t.larger[1] <- 0

  if (is.null(nparents))
    {
      if (is.null(lambda)) stop("please specify either the number of parents or the intensity of parents process")
      if (is.function(lambda)) stop("please specify the number of parents")
      if (is.numeric(lambda))
        lambda <- lambda*(areapl(s.larger)*(t.larger[2]-t.larger[1]))/(s.area*t.area)
      npar <- rpois(1,lambda)
    }
  else npar <- nparents

  if (is.null(lambda)) lambda <- (npar/(s.area*t.area))*(areapl(s.larger)*(t.larger[2]-t.larger[1]))/(s.area*t.area)
    
  if (is.null(npoints))
    {
      if (is.null(mc)) stop("please specify either the number of points to simulate or the mean number of children per parents")
      else
        {
          nchild <- matrix(rpois(npar,lambda=mc),npar,1)
          npts <- sum(nchild)
        }
    }
  else
    {
      npts <- npoints
      if (is.null(mc)) mc <- npts/npar
      nchild <- matrix(rpois(npar,lambda=mc),npar,1)
      nc <- sum(nchild)
      while (nc < npts)
        {
          nchild <- matrix(rpois(npar,lambda=mc),npar,1)
          nc <- sum(nchild)
        }
    }

  parpts <- rpp(lambda=lambda,s.region=s.larger,t.region=t.larger,npoints=npar,...)$xyt

  pattern.interm <- NULL
  ipt <- 1
  nchild2 <- rep(0,npar)
  for (ipar in 1:npar)
    {
      xpar <- parpts[ipar,1]
      ypar <- parpts[ipar,2]
      zpar <- parpts[ipar,3]
      nc <- nchild[ipar,1]
      
      if (s.distr=="uniform")
        {
          xp <- xpar+runif(nc,min=-maxrads,max=maxrads)
          yp <- ypar+runif(nc,min=-maxrads,max=maxrads)
        }
      if (s.distr=="normal")
        {
          xp <- xpar+rnorm(nc,mean=0,sd=maxrads/2)
          yp <- ypar+rnorm(nc,mean=0,sd=maxrads/2)
        }
      if (s.distr=="exponential")
        {
          xp <- xpar+rexp(nc,rate=1/maxrads)
          yp <- ypar+rexp(nc,rate=1/maxrads)
        }
      if (t.distr=="uniform")
        {
          if (infecD==TRUE) 
            zp <- zpar+runif(nc,min=-maxradt,max=maxradt)
          else
            zp <- zpar+runif(nc,min=0,max=maxradt)
        }
      if (t.distr=="normal")
        {
          if (infecD==TRUE) 
            zp <- zpar+abs(rnorm(nc,mean=0,sd=maxradt/2))
          else
            zp <- zpar+rnorm(nc,mean=0,sd=maxradt/2)
        }
      if (t.distr=="exponential")
        {
          if (infecD==TRUE) 
            zp <- zpar+rexp(nc,rate=1/maxradt)
          else
            zp <- zpar+sample(c(-1,1),1)*rexp(nc,rate=1/maxradt)
        }

      mask <- ((inout(as.points(x=xp,y=yp),poly=s.region)==T) & (zp > t.region[1] & zp < t.region[2]) & (sqrt( (((xpar-xp)^2)/maxrads^2) + (((ypar-yp)^2)/maxrads^2)) < 1))
      nchild2[ipar] <- sum(mask)
      pattern.interm <- rbind(pattern.interm,cbind(xp[mask],yp[mask],zp[mask],rep(ipar,sum(mask))))
    }

  ipt <- dim(pattern.interm)[[1]]
  if (ipt > npts)
    {
      samp <- sample(1:ipt,npts)
      pattern.interm <- pattern.interm[samp,]
    }
  while(ipt < npts)
    {
      ipar <- sample(1:npar,1)
      nchild2[ipar] <- nchild2[ipar]+1
      xpar <- parpts[ipar,1]
      ypar <- parpts[ipar,2]
      zpar <- parpts[ipar,3]
      
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
      
      if ((inout(as.points(x=xp,y=yp),poly=s.region)==T) & (zp > t.region[1] & zp < t.region[2]) & (sqrt( (((xpar-xp)^2)/maxrads^2) + (((ypar-yp)^2)/maxrads^2)) < 1))
        {
          pattern.interm <- rbind(pattern.interm,cbind(xp,yp,zp,ipar))
          ipt <- ipt+1
        }
    }
    
  pts <- pattern.interm[,1:3]
  ott<-order(pts[,3])
  pts<-pts[ott,]
  invisible(return(list(pts=pts,s.larger=s.larger,t.larger=t.larger,index.child=pattern.interm[,4],nchild=nchild2,parpts=parpts)))
}




rpcp <- function(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, mc=NULL, nsim=1, cluster="uniform", maxrad, infectious=TRUE, edge = "larger.region", ...)
{
  #
  # Simulate a space-time Poisson cluster process in a region D x T.
  # Children are simulated within a cylinder defined around the parent. 
  #
  # Requires Splancs package.
  #  
  # Arguments:
  #  
  #   s.region: two columns matrix specifying polygonal region containing
  #             all data locations. If s.region is missing, the unit square
  #             is considered.
  #
  #   t.region: vector containing the minimum and maximum values of
  #             the time interval.
  #
  #   nparents: number of parents. If missing, nparents is from a
  #             Poisson distribution with intensity lambda.
  #
  #    npoints: number of points to simulate. If NULL (default), the
  #             number of points is from a Poisson distribution with
  #             mean the double integral of lambda over s.region and
  #             t.region.
  #
  #    lambda: intensity of the parent process. Can be either a numeric
  #             value or a function. If NULL, it is constant and equal
  #             to nparents / volume of the domain.
  #
  #         mc: average number of children per parent. It is used when
  #             npoints is NULL. 
  #
  #       nsim: number of simulations to generate. Default is 1.
  #
  #    cluster: distribution of children. "uniform", "normal" and
  #             exponential are currently implemented.
  #             Either a single value if the distribution in space
  #             and time is the same, or a vector of length 2, given
  #             first the spatial distribution of children and then
  #             the temporal distribution.
  #
  #     maxrad: vector of length 2 giving the maximum spatial and temporal
  #             variation of the offspring.
  #             maxrads = maximum distance between parent and child (radius of
  #             a circle centred at the parent).
  #             maxradt = maximum time separiting parent and child.
  #             For a normal distribution of children, maxrad corresponds to 
  #             the 2 * standard deviation of location of children relative
  #             to their parent, such that children lies in the 95% IC
  #             of the normal distribution.
  #
  # infectious: If TRUE (default), corresponds to infectious diseases
  #             (times of children are always greater than parent's time).
  #
  #       edge: specify the edge correction to use, "weight", "larger.region",
  #             "without". 
  #
  #        ...: additional parameters of the intensity of the parent process.
  #
  # Value:
  #  xyt: matrix (or list of matrix if nsim>1) containing the points (x,y,t)
  #           of the simulated point process.
  #


  if (missing(cluster)) cluster <- "uniform"
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  if (missing(maxrad)) maxrad <- c(0.05,0.05)
  maxrads <- maxrad[1]
  maxradt <- maxrad[2]

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
  
  t.region <- sort(t.region)
  s.area <- areapl(s.region)
  t.area <- t.region[2]-t.region[1]
  pattern <- list()

  ni <- 1

  while(ni<=nsim)
    {
      if (edge=="larger.region")
        pattern.interm <- pcp.larger.region(s.region=s.region, t.region=t.region, nparents=nparents, npoints=npoints, lambda=lambda, mc=mc, cluster=cluster, maxrad=maxrad, infecD=infectious, ...)$pts

      if (edge=="without")
        pattern.interm <- pcp.larger.region(s.region=s.region, t.region=t.region, nparents=nparents, npoints=npoints, lambda=lambda, mc=mc, cluster=cluster, maxrad=maxrad, infecD=infectious, maxradlarger=c(0,0), ...)$pts

      if (nsim==1)
        pattern <- as.3dpoints(pattern.interm)
      else
        pattern[[ni]] <- as.3dpoints(pattern.interm)
      
      ni <- ni+1
    }
  
  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region)))
} 




