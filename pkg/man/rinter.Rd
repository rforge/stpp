\name{rinter}
\alias{rinter}
\title{Generate interaction point patterns}

\description{
Generate one (or several) realisation(s) of the inhibition or contagious process 
in a region S x T.
}

\usage{
 rinter(npoints, s.region, t.region, hs="step", gs="min", thetas=0, deltas,
        ht="step", gt="min", thetat=1, deltat, recent="all", nsim=1,
        discrete.time=FALSE, replace=FALSE, inhibition=TRUE)
}

\arguments{
  \item{npoints}{number of points to simulate. }
  \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations. 
   If \code{s.region} is missing, the unit square is considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval. 
   If \code{t.region} is missing, the interval [0,1] is considered.} 
  \item{hs, ht}{function which depends on the distance between points
  and \code{theta}. Can be chosen among "\code{step}" and "\code{gaussian}" 
  or can refer to a user defined function which only depend on d, theta, and delta
  (see details). 
  If \code{inhibition}=TRUE, \code{h} is monotone, increasing, and must tend
  to 1 when the distance tends to infinity. 0 \eqn{\leq}{<=} h(d,theta) 
  \eqn{\leq}{<=} 1. 
   Otherwise, \code{h} is monotone, decreasing, and must tend
  to 1 when the distance tends to 0.}
  \item{gs, gt}{Must be choosen among "\code{min}", "\code{max}" and 
   "\code{prod}". }
  \item{thetas, thetat}{Parameters of \code{hs} and \code{ht} functions.} 
  \item{deltas, deltat}{Spatial and temporal distance of inhibition.}
  \item{recent}{If ``\code{all}'' consider all previous events. If
  is an integer, say N, consider only  the N most recent events.}
  \item{nsim}{number of simulations to generate. Default is 1. }
  \item{discrete.time}{if TRUE, times are discrete, otherwise belong to R+.} 
  \item{replace}{Logical. If TRUE allows times repeat.}
  \item{inhibition}{Logical. If TRUE, an inhibition process is
  generated. Otherwise, it is a contagious process. }
}


#\subsubsection*{Details}
#
#Inhibition process ($k$th step):
#\begin{enumerate}
#\item Generate uniformly a location $s \in S$ and a time $t \in T$.
#\item Generate $u_s \sim {\cal U}[0,1]$ and $u_t \sim {\cal U}[0,1]$.
#\item  If $\| s -s_j \| \geq \delta_s$ for all $j=1,\dots,k-1$,
#    then set $p_s = 1$.
#
#    Otherwise compute $v_s = p_s \left( h_s
#      \left( (\| s -s_j \|)_{j=l,\dots,k-1}, \theta_s, \delta_s
#      \right), r_k \right)$
#  \item If $| t -t_j | \geq \delta_t$ for all $j=1,\dots,k-1$,
#    then set $p_t = 1$.
#
#    Otherwise compute $v_t = p_t \left( h_t \left( (| t -t_j |)_{j=l,\dots,
#          k-1}, \theta_t, \delta_t \right), r_k \right)$
#  \item If $u_s < v_s$ and $u_t < v_t$, then keep $s$ and $t$.
#\end{enumerate}
#The functions $p_s$ and $p_t$ can be choosen among ``min'', ``max''
#and ``prod'' and depend on the parameter $r_k$. This parameter allows us
#to consider either all previous events or only the most recent in the
#conditions 3 et 4.
#The functions $h_s$ and $h_t$ depend on the distance between locations
#and times respectively. They are monotone, increasing, tend
#to 1 when the distance tends to infinity and satisfy $0 \leq h(\cdot) \leq 1$.
#Currently the following functions are implemented:
#\begin{itemize}
#\item step: $h(x) = \lce
#  \begin{array}[l]{ll}
#    1,  \text{ if } x>\delta \\
#    \theta,  \text{ otherwise}
#  \end{array} \right., \ \theta \in [0,1]
#  $
#\item exponential: $h(x) = \exp(\theta x) / \max(\exp(\theta x)), \ \theta
#  \geq 0$,
#\item gaussian: $h(x) = \exp(\theta x^2) / \max(\exp(\theta x^2)), \ \theta
#  \geq 0$.
#\end{itemize}
#Note that the user can call his own functions $h_s$ and $h_t$, provided 
#that they have (spatial or temporal distance between points), 
#\verb#theta# and \verb#delta# (even if not used by the function).
#
#{\bf *** should we also allow an user defined $p$ function?}
#
#\noindent
#The simple sequential inhibition process is given using ``min'' for
#$g$ functions, ``step'' for $h$ functions and 0 for $\theta$ parameter.
#
#Contagious process ($k$th step):
#\begin{enumerate}
#\item Generate uniformly a location $s \in S$ and a time $t \in T$.
#\item Generate $u_s \sim {\cal U}[0,1]$ and $u_t \sim {\cal U}[0,1]$.
#\item If $\| s_{k-1} -s \| < \delta_s$, then set $p_s = 1$.
#  Otherwise, compute $p_s =h_s \left( \| s_{k-1} -s \|, \theta_s,
#    \delta_s \right)$.
#\item If $| t_{k-1} - t | < \delta_t$, then set $p_t=1$.
#  Otherwise, compute $p_t = h_t \left( | t_{k-1} - t|, \theta_t,
#    \delta_t \right)$.
#\item If $u_s < p_s$ and $u_t < p_t$, then keep $s$ and $t$.
#\end{enumerate}
#The functions $h_s$ and $h_t$ depend on the distance between locations
#and times respectively. They are monotone, decreasing, tend
#to 1 when the distance tends to 0 and satisfy $0 \leq h(\cdot) \leq 1$.
#$h_{cont.} (x,\theta,\delta)=h_{inhib.}(\max(x)-x,\theta,\max(x)-\delta)$.
#The simple contagious process is given using ``step'' for $h$ functions
#and 0 for $\theta$ parameters.

value{
A list containing:
\item{xyt}{matrix (or list of matrices if \code{nsim}>1)
containing the points (x,y,t) of the simulated point pattern.}
\item{s.region, t.region}{parameters passed in argument.}
}

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{animation}} and \code{\link{stani}} for plotting space-time point patterns.
 }
 
\examples{
# simple inhibition process
inh1 = rinter(npoints=200,thetas=0,deltas=0.05,thetat=0,deltat=0.001,inhibition=TRUE)
stani(inh1$xyt)

# inhibition process using hs and ht defined by the user
hs = function(d,theta,delta,mus=0.1)
{
 res=NULL
 a=(1-theta)/mus
 b=theta-a*delta
 for(i in 1:length(d))
	{	
	if (d[i]<=delta) res=c(res,theta)
	if (d[i]>(delta+mus)) res=c(res,1)
	if (d[i]>delta & d[i]<=(delta+mus)) res=c(res,a*d[i]+b)
	}
 return(res)
}
ht = function(d,theta,delta,mut=0.3)
{
 res=NULL
 a=(1-theta)/mut
 b=theta-a*delta
 for(i in 1:length(d))
	{	
	if (d[i]<=delta) res=c(res,theta)
	if (d[i]>(delta+mut)) res=c(res,1)
	if (d[i]>delta & d[i]<=(delta+mut)) res=c(res,a*d[i]+b)
	}
 return(res)
}

x=seq(0,1,length=100)
plot(x,hs(x,theta=0.2,delta=0.1,mus=0.1),xlab="",ylab="",type="l",ylim=c(0,1))
lines(x,ht(x,theta=0.1,delta=0.05,mut=0.3),col=2)
inh2 = rinter(npoints=100, hs=hs, gs="min", thetas=0.2, deltas=0.1, ht=ht, gt="min",	 thetat=0.1, deltat=0.05, inhibition=TRUE)
animation(inh2$xyt,runtime=15,cex=0.8)

# simple contagious process for given spatial and temporal regions
data(northcumbria)
cont1 = rinter(npoints=250, s.region=northcumbria, t.region=c(1,200), hs="gaussian",
 gs="min", thetas=1000, deltas=5000, ht="step", gt="min", thetat=0, deltat=10,
 recent=1, inhibition=FALSE)
animation(cont1$xyt,  s.region=cont1$s.region, t.region=cont1$t.region,
 incident="darkgreen", prevalent="green3", runtime=15, cex=0.8)
} 
