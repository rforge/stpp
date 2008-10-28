\name{stani}
\alias{stani}
\title{3D space-time data animation}
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
  \item{states}{How to display points - see Details}
}
\details{
This function requires the \code{rpanel} and \code{rgl} packages. It
uses \code{rpanel} for the sliders to control the graphics, and
\code{rgl} for its ability to do flicker-free graphics.

The sliders set the position and width of the temporal highlight
window. For 'time' slider set to time T and 'width' slider set to W,
highlighted points are those with time coordinate t such that T-W < t <
T.

How points are shown is configured with the states parameter. This is a
list of length 3 specifying how points before the time window, inside
the time window, and after the time window are displayed. Each element
is a list of parameters as would be passed to material3d() together with
a radius element. Points are drawn as spheres with the corresponding
material and radius as a fraction of the spatial span of the data.

By default the third state is invisible, and the first two states are
different. By calling with the default for states and persist=TRUE,
then the first state is set to the same as the second state. This has
the effect of showing all points at time < T with the same sphere type.

If the user specifies the states parameter, then persist is ignored. The
user can emulate the persist behaviour by specifying a states list with
identical parameters for states 1 and 2.

Note that each state element should specify all material3d parameters
used in any of the state elements. This is to make sure the parameters
are reset for each of the sets of points.

}
\value{
A list of the slider parameters when the dialog is quitted.
}
\author{Barry Rowlingson <b.rowlingson@lancaster.ac.uk>}

\keyword{ hplot }

