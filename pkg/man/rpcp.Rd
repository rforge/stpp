\name{rpcp}
\alias{rpcp}
\title{Generate Poisson cluster point patterns}

\description{
Generate one (or several) realisation(s) of the Poisson cluster process in a region 
S x T.
}

\usage{
 rpcp(s.region, t.region, nparents=NULL, npoints=NULL, lambda=NULL, 
 mc=NULL, nsim=1, cluster="uniform", maxrad, infectious=TRUE, 
 edge = "larger.region", ...)
}

\arguments{
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations.  
   If \code{s.region} is missing, the unit square is considered.} 
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing, the interval [0,1] is 
  considered.}
  \item{nparents}{number of parents. If NULL, \code{nparents} is from a
  Poisson distribution with intensity \code{lambda}.} 
  \item{npoints}{number of points to simulate. If NULL (default), the
  number of points is from a Poisson distribution with mean the double integral 
  of the intensity over \code{s.region} and \code{t.region}.} 
  \item{lambda}{intensity of the parent process. Can be either a numeric
  value, a function, or a character (see \code{\link{rpp}}). 
   If \code{NULL}, it is constant and equal to \code{nparents} / volume
  of the domain.} 
  \item{mc}{average number of children per parent. It is used when
  \code{npoints} is NULL.}  
  \item{nsim}{number of simulations to generate. }
  \item{cluster}{distribution of children: ``uniform'', ``normal'' and
  ``exponential'' are currently implemented. 
  Either a single value if the distribution in space and time is the 
   same, or a vector of length 2, giving first the spatial distribution of
  children  and then the temporal distribution.}
  \item{maxrad}{2-vector giving the maximum spatial and temporal
  variation of children. If missing, \code{maxrad}=(0.05,0.05). 
   The maximum spatial variation is the maximum distance between parent
  and children (radius of a circle centred at the parent). For a normal
  distribution of children, \code{maxrad} corresponds to the 2 * standard 
  deviation of location of children relative to their parent, such that children 
  lies in the 95 \% IC of the normal distribution. 
  The maximum temporal variation is the maximum time 
   separating parent and children.} 
  \item{infectious}{If TRUE, offspring's times are always greater than 
  parent's time).}
  \item{edge}{specify the edge correction to use "larger.region" or "without".} 
  \item{...}{additional parameters of the intensity of the parent process.} 
}

\value{
A list containing:
\item{xyt}{matrix (or list of matrices if \code{nsim}>1)
containing the points (x,y,t) of the simulated point pattern.
\code{xyt} (or any element of the list if \code{nsim}>1)) is an object 
of the class \code{stpp}.}
\item{s.region, t.region}{parameters passed in argument.}
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{animation}} and \code{\link{stani}} for plotting space-time point patterns.
 }

\examples{
# homogeneous Poisson distribution of parents

data(northcumbria)
pcp1 <- rpcp(nparents=50, npoints=500, s.region=northcumbria, t.region=c(1,365),cluster=c("normal","exponential"), maxrad=c(5000,5))
\dontrun{
animation(pcp1$xyt, s.region=pcp1$s.region, t.region=pcp1$t.region,runtime=5)
}
# inhomogeneous Poisson distribution of parents

lbda <- function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
pcp2 <- rpcp(nparents=50, npoints=500, cluster="normal", lambda=lbda, a=4000/((1-exp(-4))*(1-exp(-2))))
\dontrun{
stani(pcp2$xyt)
}

}
