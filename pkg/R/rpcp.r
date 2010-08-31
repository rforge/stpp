itexp <- function(u, m, t) { -log(1-u*(1-exp(-t*m)))/m }
rtexp <- function(n, m, t) { itexp(runif(n), m, t) }

pcp.larger.region <- function(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, mc=NULL, nsim=1, cluster="uniform", dispersion, infecD=TRUE, larger.region, tronc=1, ...)
{
  if (missing(cluster)) cluster <- "uniform"
  
  if (missing(s.region)) s.region <- matrix(c(0,0,1,1,0,1,1,0),ncol=2)
  if (missing(t.region)) t.region <- c(0,1)

  if (missing(dispersion)) dispersion <- c(0.05,0.05)
  if (length(dispersion)==1) dispersion=rep(dispersion,2)
  dispersions <- dispersion[1]
  dispersiont <- dispersion[2]

  if (missing(larger.region)) larger.region <- dispersion
  larger.regions <- larger.region[1]
  larger.regiont <- larger.region[2]
  
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

  if (larger.regions==0)
    s.larger <- s.region
  else
    {
      s.larger <- rbind(s.region+larger.regions,s.region-larger.regions,cbind(s.region[,1]+larger.regions,s.region[,2]-larger.regions),cbind(s.region[,1]-larger.regions,s.region[,2]+larger.regions))
      M <- chull(s.larger)
      s.larger <- cbind(s.larger[M,1],s.larger[M,2])
    }
  t.larger <- c(t.region[1]-larger.regiont,t.region[2]+larger.regiont)

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
          xp <- xpar+runif(nc,min=-dispersions,max=dispersions)
          yp <- ypar+runif(nc,min=-dispersions,max=dispersions)
        }
      if (s.distr=="normal")
        {
          xp <- xpar+rnorm(nc,mean=0,sd=dispersions/2)
          yp <- ypar+rnorm(nc,mean=0,sd=dispersions/2)
        }
      if (s.distr=="exponential")
        {
          xp <- xpar+rexp(nc,rate=1/dispersions)
          yp <- ypar+rexp(nc,rate=1/dispersions)
        }
      if (t.distr=="uniform")
        {
          if (infecD==TRUE) 
            zp <- zpar+runif(nc,min=-dispersiont,max=dispersiont)
          else
            zp <- zpar+runif(nc,min=0,max=dispersiont)
        }
      if (t.distr=="normal")
        {
          if (infecD==TRUE) 
            zp <- zpar+abs(rnorm(nc,mean=0,sd=dispersiont/2))
          else
            zp <- zpar+rnorm(nc,mean=0,sd=dispersiont/2)
        }
      if (t.distr=="exponential")
        {
          if (infecD==TRUE) 
            zp <- zpar+rexp(nc,rate=1/dispersiont)
          else
            zp <- zpar+sample(c(-1,1),1)*rexp(nc,rate=1/dispersiont)
        }
      if (t.distr=="trunc.exponential")
        {
          if (infecD==TRUE) 
            zp <- zpar+rtexp(nc,m=1/dispersiont,t=tronc)
          else
            zp <- zpar+sample(c(-1,1),1)*rtexp(nc,m=1/dispersiont,t=tronc)
        }

      mask <- ((inout(as.points(x=xp,y=yp),poly=s.region)==T) & (zp > t.region[1] & zp < t.region[2]) & (sqrt( (((xpar-xp)^2)/dispersions^2) + (((ypar-yp)^2)/dispersions^2)) < 1))
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
          xp <- xpar+runif(1,min=-dispersions,max=dispersions)
          yp <- ypar+runif(1,min=-dispersions,max=dispersions)
        }
      if (s.distr=="normal")
        {
          xp <- xpar+rnorm(1,mean=0,sd=dispersions/2)
          yp <- ypar+rnorm(1,mean=0,sd=dispersions/2)
        }
      if (s.distr=="exponential")
        {
          xp <- xpar+rexp(1,rate=1/dispersions)
          yp <- ypar+rexp(1,rate=1/dispersions)
        }
      if (t.distr=="uniform")
        {
          if (infecD==TRUE) 
            zp <- zpar+runif(1,min=-dispersiont,max=dispersiont)
          else
            zp <- zpar+runif(1,min=0,max=dispersiont)
        }
      if (t.distr=="normal")
        {
          if (infecD==TRUE) 
            zp <- zpar+abs(rnorm(1,mean=0,sd=dispersiont/2))
          else
            zp <- zpar+rnorm(1,mean=0,sd=dispersiont/2)
        }
      if (t.distr=="exponential")
        {
          if (infecD==TRUE) 
            zp <- zpar+rexp(1,rate=1/dispersiont)
          else
            zp <- zpar+sample(c(-1,1),1)*rexp(1,rate=1/dispersiont)
        }
      if (t.distr=="trunc.exponential")
        {
          if (infecD==TRUE) 
            zp <- zpar+rtexp(1,m=1/dispersiont,t=tronc)
          else
            zp <- zpar+sample(c(-1,1),1)*rtexp(1,m=1/dispersiont,t=tronc)
        }
      
      if ((inout(as.points(x=xp,y=yp),poly=s.region)==T) & (zp > t.region[1] & zp < t.region[2]) & (sqrt( (((xpar-xp)^2)/dispersions^2) + (((ypar-yp)^2)/dispersions^2)) < 1))
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




rpcp <- function(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, mc=NULL, nsim=1, cluster="uniform", dispersion, infectious=TRUE, edge = "larger.region", larger.region=larger.region, tronc=1,...)
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
  #     dispersion: vector of length 2 giving the maximum spatial and temporal
  #             variation of the offspring.
  #             dispersions = maximum distance between parent and child (radius of
  #             a circle centred at the parent).
  #             dispersiont = maximum time separiting parent and child.
  #             For a normal distribution of children, dispersion corresponds to 
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

  if (missing(dispersion)) dispersion <- c(0.05,0.05)
  dispersions <- dispersion[1]
  dispersiont <- dispersion[2]

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
        pattern.interm <- pcp.larger.region(s.region=s.region, t.region=t.region, nparents=nparents, npoints=npoints, lambda=lambda, mc=mc, cluster=cluster, dispersion=dispersion, infecD=infectious, larger.region=larger.region,tronc=tronc, ...)$pts

      if (edge=="without")
        pattern.interm <- pcp.larger.region(s.region=s.region, t.region=t.region, nparents=nparents, npoints=npoints, lambda=lambda, mc=mc, cluster=cluster, dispersion=dispersion, infecD=infectious, larger.region=c(0,0), tronc=tronc, ...)$pts

      if (nsim==1)
        pattern <- as.3dpoints(pattern.interm)
      else
        pattern[[ni]] <- as.3dpoints(pattern.interm)
      
      ni <- ni+1
    }
  
  invisible(return(list(xyt=pattern,s.region=s.region,t.region=t.region)))
} 




