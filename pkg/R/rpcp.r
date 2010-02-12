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




