\name{rlgcp}
\alias{rlgcp}
\title{Generate log-Gaussian Cox point patterns}

\description{
Generate one (or several) realisation(s) of the log-Gaussian cox process 
in a region S x T.
}

\usage{
 rlgcp(s.region, t.region, npoints=NULL, nsim=1, separable=TRUE, model="exponential",
       param=c(1,1,1,1,1,2), var.grf=1, mean.grf=0, replace=TRUE, nx=100, ny=100,
       nt=100, lmax=NULL, discrete.time=FALSE, exact=FALSE)
}

\arguments{
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations. If \code{s.region} is missing, the unit square is considered.} 
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing, the interval [0,1] is considered.} 
  \item{npoints}{number of points to simulate. If NULL, the
  number of points is from a Poisson distribution with mean the double integral of 	  
  lambda(s,t) over \code{s.region} and \code{t.region}. }
  \item{nsim}{number of simulations to generate. Default is 1.}
  \item{separable}{Logical. If TRUE, the covariance function of the
  Gaussian random field is separable. }
  \item{model}{vector of length 1 or 2 specifying the model(s) of
  covariance of the Gaussian random field. If \code{separable=TRUE} and 
  \code{model} is of length 2, then the elements of \code{model} define the 
  spatial and temporal covariances respectively.
   If \code{separable=TRUE} and \code{model} is of length 1, then the
  spatial and temporal covariances belongs to the same class of covariances, 
  among "exponential", "stable", "cauchy" and "wave" (see Details). 
   If \code{separable=FALSE}, \code{model} must be of length 1 and is either
  "gneiting" or "cesare" (see Details).} 
  \item{param}{vector of parameters of the covariance function (see
  Details).} 
  \item{var.grf}{variance of the Gaussian random field.}
  \item{mean.grf}{mean of the Gaussian random field.}
  \item{replace}{logical allowing times repeat.} 
  \item{nx,ny,nt}{define the size of the 3-D grid on which the intensity
  is evaluated.} 
  \item{lmax}{upper bound for the value of lambda(x,y,t).} 
  \item{discrete.time}{if TRUE, times belong to \eqn{\mathbb{N}}{N},
  otherwise belong to \eqn{\mathbb{R}^+}{R^+}.} 
}

\details{
We implemented stationary, isotropic spatio-temporal covariance functions.

\emph{Separable covariance functions}

\deqn{c(h,t) = c_s(\| h \|) \, c_t(|t|) , h \in S, t \in T}{c(h,t)
=c_s(||h||) c_t(|t|), h in S, t in T}

The purely spatial and purely temporal covariance functions can be:
\itemize{
\item Exponential: \eqn{c(r) = \exp(-r)}{c(r)=exp(-r)},
\item Stable: \eqn{c(r) = \exp(-r^\alpha)}{c(r)=exp(-r \alpha)}, 
\eqn{\alpha \in [0,2]}{\alpha \in [0,2]},
\item Cauchy: \eqn{c(r) = (1+r^2)^{-\alpha}}{c(r)=(1+r^2)^(-\alpha)},
\eqn{\alpha >0}{\alpha > 0},
\item Wave: \eqn{c(r) = \frac{\sin(r)}{r}}{c(r)=sin(r)/r} if \eqn{r>0}{r>0}, 
 \eqn{c(0)=1}{c(0)=1},
\item Matèrn: \eqn{c(r) = \frac{(r/\beta)^\alpha}{2^{\alpha-1}\Gamma(\alpha)}
  {\cal K}_{\alpha}(r/\beta)}{c(r)={(r/\beta)^\alpha}/{2^(\alpha-1) \Gamma(\alpha)}
  K_\alpha(r/\beta)}, \eqn{\alpha > 0}{\alpha>0} and \eqn{\beta > 0}{\beta>0}.

\eqn{{\cal K}_{\alpha}}{K_\alpha} is the modified Bessel function of second kind:
\deqn{{\cal K}_{\alpha}(x) = \frac{\pi}{2} \frac{I_{-\alpha}(x) -
    I_{\alpha}(x)}{\sin(\pi \alpha)},}{K_\alpha(x) = (\pi/2) * (I_(-\alpha)(x) 
    - I_\alpha(x))/sin(\pi \alpha),} with
  \eqn{I_{\alpha}(x) = \left( \frac{x}{2} \right)^{\alpha} \sum_{k=0}^\infty
  \frac{1}{k! \Gamma(\alpha+k+1)} \left( \frac{x}{2} \right)^{2k}}{I_\alpha(x) 
  = (x/2)^\alpha sum_{k=0}^\infty 1/(k! \Gamma(\alpha+k+1)) (x/2)^(2k)}
The parameters \eqn{alpha_1}{\alpha_1} and \eqn{\alpha_2}{\alpha_2} correspond 
to the parameters of the spatial and temporal covariance respectively.
}


\emph{Non-separable covariance functions}

The spatio-temporal covariance function can be:
\itemize{
\item gneiting: \eqn{c(h,t) = \psi (t/\beta_2)^{-\alpha_6} \phi \left(\frac{h/
      \beta_1}{\psi (t/\beta_2)} \right)}{c(h,t)=\psi(t/\beta_2)^(-\alpha_6) 
      \phi( (h/\beta_1)/\psi(t/\beta_2) )}, \eqn{\beta_1, \beta_2 >0}{\beta_1, 
      \beta_2 >0},
      \itemize{
  \item If \eqn{\alpha_2=1}{\alpha_2=1}, \eqn{\phi(r)}{\phi(r)} is the Stable model.
  \item if \eqn{\alpha_2=2}{\alpha_2=2}, \eqn{\phi(r)}{\phi(r)} is the Cauchy model.
  \item If \eqn{\alpha_2=3}{\alpha_2=3}, \eqn{\phi(r)}{\phi(r)} is the Matèrn model.
  \item If \eqn{\alpha_5=1}{\alpha_5=1}, \eqn{\psi^2(r) = (r^{\alpha_3}+
  1)^{\alpha_4}}{\psi^2(r) = (r^\alpha_3+1)^\alpha_4},
  \item If \eqn{\alpha_5=2}{\alpha_5=2}, \eqn{\psi^2(r) = (\alpha_4^{-1} r^{\alpha_3}
    +1)/(r^{\alpha_3}+1)}{\psi^2(r) = (\alpha_4^(-1) r^\alpha_3 +1)/(r^\alpha_3+1)},
  \item If \eqn{\alpha_5=3}{\alpha_5=3}, \eqn{\psi^2(r) = -\log(r^{\alpha_3} +
   1/{\alpha_4})/\log {\alpha_4}}{\psi^2(r) = -log(r^\alpha_3+1/\alpha_4)/
   log(\alpha_4)},
                }
        
  The parameter \eqn{\alpha_1}{\alpha_1} is the respective parameter for the model of
  \eqn{\phi(\cdot)}{\phi(.)}, \eqn{\alpha_3 \in (0,2]}{\alpha_3 in (0,2]}, 
  \eqn{\alpha_4 \in (0,1]}{\alpha_4 in (0,1]} and \eqn{\alpha_6 \geq 2}{\alpha_6 >= 2}.
\item cesare: \eqn{c(h,t) = \left( 1 + (h/\beta_1)^{\alpha_1} +
    (t/\beta_2)^{\alpha_2} \right)^{-\alpha_3}}{c(h,t) = (1 + (h/\beta_1)^\alpha_1 +
    (t/\beta_2)^\alpha_2)^(-\alpha_3)}, \eqn{\beta_1, \beta_2 >0}{\beta_1, 
    \beta_2 >0}, \eqn{\alpha_1, \alpha_2 \in [1,2]}{\alpha_1, \alpha_2 in [1,2]}
     and \eqn{\alpha_3 \geq 3/2}{\alpha_3 >= 3/2}.
	}
}

\value{
A list containing:
\item{xyt}{matrix (or list of matrices if \code{nsim}>1)
containing the points (x,y,t) of the simulated point pattern.}
\item{s.region, t.region}{parameters passed in argument.}
\item{Lambda}{nx * ny * nt array (or list of array if \code{nsim}>1) 
of the intensity.}
}

\references{
Chan, G. and Wood A. (1997).  An algorithm for simulating stationary
Gaussian random fields. Applied Statistics, Algorithm Section, 46, 171--181.

Gneiting T. (2002). Nonseparable, stationary covariance functions
for space-time data. Journal of the American Statistical Association,
97, 590--600.
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{animation}} and \code{\link{stani}} for plotting space-time point patterns.
 }

\examples{
# non separable covariance function: 
lgcp1 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=FALSE, model="gneiting",
 param=c(1,1,1,1,1,2), var.grf=1, mean.grf=0)
N <- lgcp1$Lambda[,,1];for(j in 2:(dim(lgcp1$Lambda)[3])){N <- N+lgcp1$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp1$xyt, cex=0.8, runtime=10, add=TRUE, prevalent="orange")

# separable covariance function: 
lgcp2 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=TRUE, model="exponential",
 param=c(1,1,1,1,1,2), var.grf=2, mean.grf=-0.5*2)
N <- lgcp2$Lambda[,,1];for(j in 2:(dim(lgcp2$Lambda)[3])){N <- N+lgcp2$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp2$xyt, cex=0.8, pch=20, runtime=10, add=TRUE, prevalent="orange")
}
