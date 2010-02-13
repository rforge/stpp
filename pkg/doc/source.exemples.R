

library(stpp)

data(fmd)
data(northcumbria) 

fmd=as.3dpoints(fmd)
plot(fmd,s.region=northcumbria)

animation(fmd,runtime=10,cex=0.5,s.region=northcumbria)

library(rgl)
library(rpanel)

stani(fmd,bgpoly=northcumbria,bgframe=FALSE,twid=1) 
# try with the first 100 points if too long :
# stani(fmd[1:100,],bgpoly=northcumbria,bgframe=FALSE,twid=1) 


###############
#
#	hpp
#
###############

hpp1 = rpp(lambda=200, nsim=5, replace=FALSE)
stani(hpp1$xyt[[2]])

data(northcumbria)
hpp2 = rpp(npoints=1000, s.region=northcumbria, t.region=c(1,500), discrete.time=TRUE)
polymap(northcumbria)
animation(hpp2$xyt, s.region=hpp2$s.region, add=TRUE)


###############
#
#	ipp
#
###############

lbda1 = function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
ipp1 = rpp(lambda=lbda1, npoints=200, a=1600/((1-exp(-4))*(1-exp(-2))))
stani(ipp1$xyt)

data(fmd)
data(northcumbria)
h = mse2d(as.points(fmd[,1:2]), northcumbria, nsmse=30, range=3000)
h = h$h[which.min(h$mse)]
Ls = kernel2d(as.points(fmd[,1:2]), northcumbria, h, nx=100, ny=100)
Lt = dim(fmd)[1]*density(fmd[,3], n=200)$y
ipp2 = rpp(lambda="m", Lambda=Ls$z, mut=Lt, s.region=northcumbria, 
  t.region=c(1,200), discrete.time=TRUE)
image(Ls$x, Ls$y, Ls$z, col=grey((1000:1)/1000)); polygon(northcumbria)
animation(ipp2$xyt, add=TRUE, cex=0.5, runtime=15)


###############
#
#	pcp
#
###############

data(northcumbria)
pcp1 <- rpcp(nparents=50, mc=10, s.region=northcumbria, t.region=c(1,365),
  cluster=c("normal","exponential"), maxrad=c(5000,5))
animation(pcp1$xyt, s.region=pcp1$s.region, t.region=pcp1$t.region,runtime=5)

lbda <- function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
pcp2 <- rpcp(nparents=50, npoints=250, cluster="normal", lambda=lbda, 
  a=2000/((1-exp(-4))*(1-exp(-2))))
stani(pcp2$xyt)

# inhibition

inh1 = rinter(npoints=200, thetas=0, deltas=0.05, thetat=0, deltat=0.001, inhibition=TRUE)
stani(inh1$xyt)

######################
#
#	contagious
#
######################

data(northcumbria)
cont1 = rinter(npoints=250, s.region=northcumbria, t.region=c(1,200), 
  thetas=1000, deltas=5000, thetat=0, deltat=10, recent=1, inhibition=FALSE)
animation(cont1$xyt,  s.region=cont1$s.region, t.region=cont1$t.region,
  incident="red", prevalent="lightgreen", runtime=15, cex=0.8)

# defining and plotting hs and ht
hs = function(d,theta,delta,mus=0.1){
 res=NULL
 a=(1-theta)/mus
 b=theta-a*delta
 for(i in 1:length(d))
	{	
	if (d[i]<=delta) res=c(res,theta)
	if (d[i]>(delta+mus)) res=c(res,1)
	if (d[i]>delta & d[i]<=(delta+mus)) res=c(res,a*d[i]+b)
	}
 return(res)}

ht = function(d,theta,delta,mut=0.3){
 res=NULL
 a=(1-theta)/mut
 b=theta-a*delta
 for(i in 1:length(d))
	{	
	if (d[i]<=delta) res=c(res,theta)
	if (d[i]>(delta+mut)) res=c(res,1)
	if (d[i]>delta & d[i]<=(delta+mut)) res=c(res,a*d[i]+b)
	}
 return(res)}

d=seq(0,1,length=100)
plot(d, hs(d,0.2,0.1,0.1), xlab="", ylab="", type="l", ylim=c(0,1), lwd=2, las=1)
lines(d, ht(d,0.1,0.05,0.3), col=2, lwd=2)
legend("bottomright", col=1:2, lty=1, lwd=2, bty="n", cex=2, 
  legend=c(expression(h[s]),expression(h[t])))

# generating the inhibition process
inh2 = rinter(npoints=100, hs=hs, gs="min", thetas=0.2, deltas=0.1, ht=ht, 
	gt="min", thetat=0.1, deltat=0.05, inhibition=TRUE)
animation(inh2$xyt,runtime=15,cex=0.8)


######################
#
#	infectious
#
#######################

inf1 = rinfec(npoints=100, alpha=0.1, beta=0.6, gamma=0.5, maxrad=c(0.075,0.5),
  t.region=c(0,50), s.distr="uniform", t.distr="uniform", h="step", g="min",
  recent="all", inhibition=TRUE)
animation(inf1$xyt, cex=0.8, runtime=10)

data(fmd)
data(northcumbria)
h = mse2d(as.points(fmd[,1:2]), northcumbria, nsmse=30, range=3000)
h = h$h[which.min(h$mse)]
Ls = kernel2d(as.points(fmd[,1:2]), northcumbria, h, nx=50, ny=50)
inf2 = rinfec(npoints=100, alpha=4, beta=0.6, gamma=20, maxrad=c(12000,20), 
  s.region=northcumbria, t.region=c(1,2000), s.distr="poisson", t.distr="uniform", 
  h="step", g="min", recent=1, lambda=Ls$z, inhibition=FALSE)
image(Ls$x, Ls$y, Ls$z, col=grey((1000:1)/1000)); polygon(northcumbria,lwd=2)
animation(inf2$xyt, add=TRUE, cex=0.7, runtime=15, prevalent=4)



###############
#
#	lgcp
#
###############

lgcp1 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=FALSE, 
  model="gneiting", param=c(1,1,1,1,1,2), var.grf=1, mean.grf=0)
N <- lgcp1$Lambda[,,1]
for(j in 2:(dim(lgcp1$Lambda)[3])){N <- N+lgcp1$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp1$xyt, cex=0.8, runtime=10, add=TRUE, prevalent="orange")

lgcp2 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=TRUE,
 model="exponential",param=c(1,1,1,1,1,2), var.grf=2, mean.grf=-0.5*2)
N <- lgcp2$Lambda[,,1]
for(j in 2:(dim(lgcp2$Lambda)[3])){N <- N+lgcp2$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp2$xyt, cex=0.8, runtime=10, add=TRUE, prevalent="orange")
