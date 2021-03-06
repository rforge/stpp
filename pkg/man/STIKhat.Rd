\name{STIKhat} \alias{STIKhat}
\title{Estimation of the Space-Time Inhomogeneous K-function}

\description{Compute an estimate of the Space-Time
Inhomogeneous K-function.}

\usage{STIKhat(xyt, s.region, t.region, dist, times, lambda,
correction=TRUE, infectious=TRUE) }

\arguments{
  \item{xyt}{coordinates and times (x,y,t) of the point pattern.}
    \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations.
   If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval.
   If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  \item{dist}{vector of distances u at which \eqn{K_{ST}(u,v)}{K(u,v)} is computed.}
   \item{times}{vector of times v at which \eqn{K_{ST}(u,v)}{K(u,v)} is computed.}
    \item{lambda}{vector of values of the space-time intensity function
  evaluated at the points (x,y,t) in SxT.
  If \code{lambda} is missing, the estimate of the space-time
K-function is computed as for the homogeneous case (Diggle et
al., 1995), i.e. considering \eqn{n/|S \times T|}{n/|SxT|} as
an estimate of the space-time intensity.}
\item{correction}{logical value. If \code{TRUE}, spatial (Ripley's) and temporal
edge corrections are used.}
\item{infectious}{logical value. If \code{TRUE}, only future events
  are considered. See Details.}
 }

\details{Gabriel and Diggle (2009) propose the following
approximately unbiased estimator for the STIK-function, based
on data giving the locations of events \eqn{x_i: i=1,...n}{xi:
i=1...,n} on a spatio-temporal region \eqn{S \times
T}{SxT},where S is an arbitrary polygon and
\eqn{T=[T_0,T_1]}{T=[T0,T1]}: \deqn{
  \widehat K_{ST}(u,v) =  \frac{1}{|S \times T|} \frac{n}{n_v} \sum_{i=1}^{n_v}
  \sum_{j=1; j > i}^{n_v} \frac{1}{w_{ij}} \frac{1}{\lambda(x_i) \lambda(x_j)}{\bf 1}_{\left\lbrace u_{ij} \leq u
  \right\rbrace} {\bf 1}_{\left\lbrace t_j - t_i \leq v \right\rbrace}}{K(u,v) = 1/|SxT| n/nv sum_{i=1,...,nv}
  sum_{j=1,...,nv; j > i} 1/wij 1/(lambda(x_i)lambda(x_j)) 1{uij <=
  u} 1{tj - ti <= v}}
In this equation, \eqn{\lambda(x_i)}{lambda(xi)} is the
intensity at \eqn{x_i = (s_i,t_i)}{xi=(si,ti)} and the
\eqn{x_i}{xi} are ordered so that \eqn{t_i < t_{i+1}}{ti <
t(i+1)}, with ties due to round-off error broken by randomly
unrounding if necessary. To deal with temporal edge-effects,
for each v, \eqn{n_v}{nv} denotes the number of events for
which \eqn{t_i \leq T_1 -v}{ti <= T1-v}. To deal with spatial
edge-effects, we use Ripley's method, in which
\eqn{w_{ij}}{wij} is the proportion of the circle centered on
\eqn{s_i}{si} and passing through \eqn{s_j}{sj}, i.e. of radius
\eqn{u_{ij}= \| s_i - s_j\|}{uij=||si-sj||}, that lies inside
S.

If \code{lambda} is missing in argument, \code{STIKhat}
computes an estimate of the space-time (homogeneous)
K-function: \deqn{\widehat{K}'_{ST}(u,v) = \frac{|S \times
T|}{n_v(n-1)} \sum_{i=1}^{n_v}
  \sum_{j=1; j > i}^{n_v} \frac{1}{w_{ij}}{\bf 1}_{\left\lbrace u_{ij} \leq u \right\rbrace} {\bf 1}_{\left\lbrace t_j - t_i \leq v
  \right\rbrace}}{K'(u,v) = |SxT|/(nv(n-1)) sum_{i=1,...,nv} sum_{j=1,...,nv; j>i} 1/wij 1{uij <=
  u} 1{tj - ti <= v}}

If parameter \code{infectious = FALSE}, both future and past
events are considered and the estimator is: \deqn{
  \widehat K^*_{ST}(u,v) =  \frac{1}{|S \times T|} \sum_{i=1}^{n}
  \sum_{j \neq j}^{n} \frac{1}{w_{ij}} \frac{1}{\lambda(x_i) \lambda(x_j)}{\bf 1}_{\left\lbrace u_{ij} \leq u
  \right\rbrace} {\bf 1}_{\left\lbrace |t_j - t_i| \leq v \right\rbrace}}{K^*(u,v) = 1/|SxT| sum_{i=1,...,n}
  sum_{j=1,...,n; j \neq j} 1/wij 1/(lambda(xi)lambda(xj)) 1{uij <=
  u} 1{|tj - ti| <= v}}
  where \eqn{v_{ij}}{vij} is equal to 1 if both ends of the
  interval of length \eqn{2 |t_i - t_j|}{2|ti-tj|} centred at \eqn{t_i}{ti} lie within
  T and 1/2 otherwise.
}

 \value{ A list containing:

 \item{Khat}{ndist x ntimes matrix containing values of \eqn{\hat
 K_{ST}(u,v)}{K(u,v)}.}
\item{Ktheo}{ndist x ntimes matrix containing theoretical
values for a Poisson process; \eqn{\pi u^2 v}{pi u^2 v} for K
and \eqn{2 \pi u^2 v}{2 pi u^2 v}) for K^*.}

 \item{dist, times, infectious}{parameters passed in argument.}
 }

\references{ Gabriel E., Diggle P. (2009) Second-order analysis
of inhomogeneous spatio-temporal point process data. Statistica
Neerlandica, 63, 43--51.

Baddeley A., Moller J. and Waagepetersen R. (2000). Non- and
semi-parametric estimation of interaction in inhomogeneous
point patterns. Statistica Neerlandica, 54, 329--350.

Diggle P. , Chedwynd A., Haggkvist R. and Morris S. (1995).
Second-order analysis of space-time clustering. Statistical
Methods in Medical Research, 4, 124--136. }

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> }


\examples{\dontrun{
data(fmd)
data(northcumbria)
FMD<-as.3dpoints(fmd[,1]/1000,fmd[,2]/1000,fmd[,3])
Northcumbria=northcumbria/1000

# estimation of the temporal intensity
Mt<-density(FMD[,3],n=1000)
mut<-Mt$y[findInterval(FMD[,3],Mt$x)]*dim(FMD)[1]

# estimation of the spatial intensity
h<-mse2d(as.points(FMD[,1:2]), Northcumbria, nsmse=50, range=4)
h<-h$h[which.min(h$mse)]
Ms<-kernel2d(as.points(FMD[,1:2]), Northcumbria, h, nx=5000, ny=5000)
atx<-findInterval(x=FMD[,1],vec=Ms$x)
aty<-findInterval(x=FMD[,2],vec=Ms$y)
mhat<-NULL
for(i in 1:length(atx)) mhat<-c(mhat,Ms$z[atx[i],aty[i]])

# estimation of the STIK function
u <- seq(0,10,by=1)
v <- seq(0,15,by=1)
stik <- STIKhat(xyt=FMD, s.region=northcumbria/1000,t.region=c(1,200), 
lambda=mhat*mut/dim(FMD)[1], dist=u, times=v, infectious=TRUE)

# plotting the estimation
plotK(stik)
plotK(stik,persp=T,theta=-65,phi=35)
}}
