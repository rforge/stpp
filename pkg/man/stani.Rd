\name{stani}
\alias{stani}
\title{Space-time data animation}
\description{
 Displays (x,y,t) point data and enables dynamic highlighting of time slices.
}
\usage{
stani(xyt, tlim = range(xyt[, 3], na.rm = TRUE), twid = diff(tlim)/20, persist = FALSE, states)
}
\arguments{
  \item{xyt}{A 3-column matrix of x,y,t coordinates}
  \item{tlim}{A two-element vector of upper and lower time limits}
  \item{twid}{The initial time window width}
  \item{persist}{Whether to display points before time window}
  \item{states}{How to display points}
}
\details{
This function requires the \code{rpanel} and \code{rgl} packages. It
uses \code{rpanel} for the sliders to control the graphics, and
\code{rgl} for its ability to do flicker-free graphics.



}
\value{
None
}
\author{Barry Rowlingson <b.rowlingson@lancaster.ac.uk>}

\keyword{ hplot }

