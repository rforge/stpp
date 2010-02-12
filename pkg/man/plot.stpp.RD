\name{plot.stpp}
\alias{plot.stpp}
\title{Plot for spatio-temporal point objects}

\description{
This function plot xy-locations and cumulative distribution of the times.}

\usage{
## S3 method for class 'stpp':
plot(p, s.region=NULL, t.region=NULL)
}

\arguments{
  \item{p}{any object of class \code{stpp} in spatio-temporal point format.}
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations. If \code{s.region} is missing, the default limits are considered.} 
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing,  the default limits are considered.}
 }

\value{
None

\seealso{
 \code{\link{as.3dpoints}} for creating data in spatio-temporal point format.
 }

\author{Edith Gabriel <edith.gabriel@univ-avignon.fr>.}
