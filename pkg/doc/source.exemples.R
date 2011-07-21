
library(stpp)
data(fmd)
data(northcumbria) 

fmd=as.3dpoints(fmd)
source("stan2.R")

xyt=cbind(runif(100),runif(100),runif(100))
stan(xyt)



library(stpp)

data(fmd)
data(northcumbria) 

fmd=as.3dpoints(fmd)
plot(fmd,s.region=northcumbria,pch=19,cex=0.5)

plot(fmd,s.region=northcumbria,pch=19,mark=TRUE,mark.cexmax=0.8)


animation(fmd,runtime=10,cex=0.5,s.region=northcumbria)

stan(fmd,bgpoly=northcumbria,axes=FALSE) 

##########################
##				##
## 	STIK function	##
##				##
##########################
data(fmd)
data(northcumbria)
FMD=as.3dpoints(fmd[,1]/1000,fmd[,2]/1000,fmd[,3])

# estimation of the temporal intensity 
M=density(FMD[,3],kernel="gaussian",n=200)
mut=M$y[FMD[,3]]*dim(fmd)[1]

# estimation of the spatial intensity
h = mse2d(pts=FMD[,1:2], poly=northcumbria/1000, nsmse=100, range=5)
hs=h$h[h$mse==min(h$mse)]
require(spatialkernel)
mhat <- lambdahat(pts=as.points(FMD[,1:2]), h=hs, gpts=as.points(FMD[,1:2]), poly = northcumbria/1000, edge = TRUE)$lambda

# estimation of the STIK function
u <- seq(0,10,by=1)
v <- seq(0,15,by=1)
stik <- STIKhat(xyt=FMD, s.region=northcumbria/1000,t.region=c(1,200), lambda=mhat*mut/dim(fmd)[1], dist=u, times=v, infectious=T)

# plotting the estimation
plotK(stik)
plotK(stik,persp=T,theta=-65,phi=35)



FMD=fmd
FMD[,1:2]=fmd[,1:2]/1000

M=density(FMD[,3],kernel="gaussian",n=200)
mut=M$y[FMD[,3]]*648

hist(FMD[,3],probability=T)
lines(M$y,col=3)

h = mse2d(pts=FMD[,1:2], poly=northcumbria/1000, nsmse=100, range=5)
hs=h$h[h$mse==min(h$mse)]
library(spatialkernel)
setkernel("gaussian")
mhat <- lambdahat(pts=as.points(FMD[,1:2]), h=hs, gpts=as.points(FMD[,1:2]), poly = northcumbria/1000, edge = TRUE)$lambda
mhat2 <- lambdahat(pts=as.points(FMD[,1:2]), h=hs2, gpts=as.points(FMD[,1:2]), poly = northcumbria/1000, edge = TRUE)$lambda


u <- seq(0,10,by=1)
v <- seq(0,15,by=1)
stik <- STIKhat(xyt=FMD, s.region=northcumbria/1000,t.region=c(1,200), lambda=mhat*mut/dim(fmd)[1], dist=u, times=v, infectious=T)
stik2 <- STIKhat(xyt=FMD, s.region=northcumbria/1000,t.region=c(1,200), lambda=mhat2*mut/dim(fmd)[1], dist=u, times=v, infectious=T)


contour(u, v, sqrt(stik$Khat-stik$Ktheo), axes=F,xlab="",ylab="",xlim=c(0.5,max(u)),ylim=c(0.5,max(v)))
box(lwd=2)
at <- axTicks(1)
axis(1,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)])
axis(1,at=at[length(at)],labels="u")
at <- axTicks(2)
axis(2,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)],las=1)
axis(2,at=at[length(at)],labels="v",las=1)

persp(u,v,sqrt(stik$Khat-stik$Ktheo),theta=-45,phi=50)


# Test for spatio-temporal clustering
# ***********************************

# estimation of the spatio-temporal intensity 

M=density(FMD[,3],kernel="gaussian",n=200)
LT=M$y*dim(FMD)[1]

h = mse2d(pts=FMD[,1:2], poly=northcumbria/1000, nsmse=50, range=5)
hs=h$h[h$mse==min(h$mse)]
xygrid <- make.grid(100,100,northcumbria/1000)
require(spatialkernel)
setkernel("gaussian")
M <- lambdahat(pts=as.points(FMD[,1:2]), h=hs, gpts=cbind(as.vector(xygrid$X),as.vector(xygrid$Y)), poly = northcumbria/1000, edge = TRUE)$lambda
LS <- matrix(M,ncol=100,byrow=T)
LS[xygrid$mask==F] <- NaN
image(xygrid$x,xygrid$y,LS,col=grey((1000:1)/1000),xlab="",ylab="");polygon(northcumbria/1000)

M <- lambdahat(pts=as.points(FMD[,1:2]), h=hs2, gpts=cbind(as.vector(xygrid$X),as.vector(xygrid$Y)), poly = northcumbria/1000, edge = TRUE)$lambda
LS2 <- matrix(M,ncol=100,byrow=T)
LS2[xygrid$mask==F] <- NaN
image(xygrid$x,xygrid$y,LS2,col=grey((1000:1)/1000),xlab="",ylab="");polygon(northcumbria/1000)

Lst=array(0,dim=c(100,100,200))
for(k in 1:200) Lst[,,k] <- LS*LT[k]/dim(fmd)[1]

Lst2=array(0,dim=c(100,100,200))
for(k in 1:200) Lst2[,,k] <- LS2*LT[k]/dim(fmd)[1]


# Simulation of IPP with intensity LSxLT

simH0 <- rpp(lambda="m", Lambda=Lst, s.region=northcumbria/1000, t.region=c(1,200), discrete.time=TRUE, nsim=39)
simH02 <- rpp(lambda="m", Lambda=Lst2, s.region=northcumbria/1000, t.region=c(1,200), discrete.time=TRUE, nsim=39)

sim=simH02

# Envelopes

nsim <- length(simH0$xyt)
stikH02 <- array(0,dim=c(length(u),length(v),nsim))
for(i in 1:nsim)
  {
    muhat0 <- LT[sim$xyt[[i]][,3]]
    mhat0 <- lambdahat(pts=as.points(FMD[,1:2]), h=hs, gpts=as.points(sim$xyt[[i]][,1:2]), poly = northcumbria/1000, edge = TRUE)$lambda	
    stikH02[,,i] <- STIKhat(xyt=sim$xyt[[i]], s.region=northcumbria/1000,t.region=c(1,200), lambda=mhat0*muhat0/dim(simH0$xyt[[i]])[1], dist=u, times=v)$Khat
}

Kenv2 <- array(0,dim=c(length(u),length(v),2))
E0 <- matrix(0,nrow=length(u),ncol=length(v))
V0 <- matrix(0,nrow=length(u),ncol=length(v))
for(i in 1:length(u))
  {
    for(j in 1:length(v))
      {
        E0 <- mean(stikH0[i,j,])
        V0 <- var(stikH0[i,j,])
        M <- quantile(stikH0[i,j,],c(0.025,0.975))
        Kenv2[i,j,] <- M
      }
  }

# p-value

TclustH0 <- rep(0,nsim)
for(i in 1:nsim)
  {
     TclustH0[i] <- sum((stikH0[,,i]-E0)/sqrt(V0),na.rm=T)
  }
Tdata <- sum((stik$Khat-E0)/sqrt(V0),na.rm=T)
pclust <- (sum(Tdata < TclustH0) + 1)/(nsim + 1)

KH0=stik$Ktheo
par(mar=c(2,2,2,0),mgp=c(3,1,0),lwd=2,font=1)
mask2 <- matrix(0,ncol=length(v),nrow=length(u))
for(i in 1:length(u))
  {
    for(j in 1:length(v))
      {
        if (stik2$Khat[i,j] > Kenv2[i,j,2]) 
          mask2[i,j] <- 1
        if (stik2$Khat[i,j] < Kenv2[i,j,1])
          mask2[i,j] <- -1
      }
  }

image(u, v, mask, axes=F,col=c("grey","white","black"),xlab="",ylab="",zlim=c(-1,1),xlim=c(0.5,max(u)),ylim=c(0.5,max(v)))
box(lwd=2)
at <- axTicks(1)
axis(1,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)])
axis(1,at=at[length(at)],labels="u")
at <- axTicks(2)
axis(2,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)],las=1)
axis(2,at=at[length(at)],labels="v",las=1)


###############
#
#	hpp
#
###############

hpp1 = rpp(lambda=200, nsim=5, replace=FALSE)
plot(hpp1$xyt[[2]],pch=19)
animation(hpp1$xyt[[2]])
stan(hpp1$xyt[[2]])

data(northcumbria)
hpp2 = rpp(npoints=1000, s.region=northcumbria, t.region=c(1,500), discrete.time=TRUE)
plot(hpp2$xyt,pch=19,s.region=hpp2$s.region,cex=0.25)
par(mfrow=c(1,1))
polymap(northcumbria)
animation(hpp2$xyt, s.region=hpp2$s.region, add=TRUE)


###############
#
#	ipp
#
###############

lbda1 = function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
ipp1 = rpp(lambda=lbda1, npoints=200, a=1600/((1-exp(-4))*(1-exp(-2))))
plot(ipp1$xyt,pch=19)
animation(ipp1$xyt)
stan(ipp1$xyt)

data(fmd)
data(northcumbria)
h = mse2d(as.points(fmd[,1:2]), northcumbria, nsmse=30, range=3000)
h = h$h[which.min(h$mse)]
Ls = kernel2d(as.points(fmd[,1:2]), northcumbria, h, nx=100, ny=100)
Lt = dim(fmd)[1]*density(fmd[,3], n=200)$y
Lst=array(0,dim=c(100,100,200))
for(k in 1:200) Lst[,,k] <- Ls$z*Lt[k]/dim(fmd)[1]
ipp2 = rpp(lambda="m", Lambda=Lst, s.region=northcumbria, 
  t.region=c(1,200), discrete.time=TRUE)
            
plot(ipp2$xyt,pch=19,s.region=ipp2$s.region,cex=0.6)
par(mfrow=c(1,1))
image(Ls$x, Ls$y, Ls$z, col=grey((1000:1)/1000)); polygon(northcumbria)
animation(ipp2$xyt, add=TRUE, cex=0.5, runtime=15)


###############
#
#	pcp
#
###############

data(northcumbria)
pcp1 <- rpcp(nparents=50, mc=10, s.region=northcumbria, t.region=c(1,365),
  cluster=c("normal","exponential"), dispersion=c(5000,5))
plot(pcp1$xyt,pch=19,s.region=pcp1$s.region,cex=0.3)
animation(pcp1$xyt, s.region=pcp1$s.region, t.region=pcp1$t.region,runtime=5)





lbda <- function(x,y,t,a){a*exp(-4*y) * exp(-2*t)}
pcp2 <- rpcp(nparents=50, npoints=250, cluster="normal", lambda=lbda, 
  a=2000/((1-exp(-4))*(1-exp(-2))))
plot(pcp2$xyt,pch=19)
animation(pcp2$xyt,runtime=5)
stan(pcp2$xyt)




######################
#
#	contagious
#
######################

data(northcumbria)
cont1 = rinter(npoints=250, s.region=northcumbria, t.region=c(1,200), thetas=0,
 deltas=7500, thetat=0, deltat=10, recent=1, inhibition=FALSE)
plot(cont1$xyt,pch=19,s.region=cont1$s.region,mark=TRUE,mark.col=4)
animation(cont1$xyt,  s.region=cont1$s.region, t.region=cont1$t.region,
  incident="red", prevalent="lightgreen", runtime=15, cex=0.8)

# inhibition

inh1 = rinter(npoints=200, thetas=0, deltas=0.05, thetat=0, deltat=0.001, inhibition=TRUE)
plot(inh1$xyt,pch=19)
animation(inh1$xyt,runtime=5)
stan(inh1$xyt)

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
plot(inh2$xyt,pch=19)
animation(inh2$xyt,runtime=15,cex=0.8)







######################
#
#	infectious
#
#######################

inf1 = rinfec(npoints=100, alpha=0.1, beta=0.6, gamma=0.5, maxrad=c(0.075,0.5),
  t.region=c(0,50), s.distr="uniform", t.distr="uniform", h="step", g="min",
  recent="all", inhibition=TRUE)
plot(inf1$xyt,pch=19)
animation(inf1$xyt, cex=0.8, runtime=10)

data(fmd)
data(northcumbria)
h = mse2d(as.points(fmd[,1:2]), northcumbria, nsmse=30, range=3000)
h = h$h[which.min(h$mse)]
Ls = kernel2d(as.points(fmd[,1:2]), northcumbria, h, nx=50, ny=50)
inf2 = rinfec(npoints=200, alpha=4, beta=0.6, gamma=20, maxrad=c(12000,20), 
  s.region=northcumbria, t.region=c(1,2000), s.distr="poisson", t.distr="uniform", 
  h="step", g="min", recent=1, lambda=Ls$z, inhibition=FALSE)
plot(inf2$xyt,pch=19,s.region=inf2$s.region,mark=TRUE)
par(mfrow=c(1,1))
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
plot(lgcp1$xyt,pch=19)
par(mfrow=c(1,1))
image(N,col=grey((1000:1)/1000));box()
animation(lgcp1$xyt, cex=0.8, runtime=10, add=TRUE, prevalent="orange")

lgcp2 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=TRUE,
 model="exponential",param=c(1,1,1,1,1,2), var.grf=2, mean.grf=-0.5*2)
plot(lgcp2$xyt,pch=19)
par(mfrow=c(1,1))
N <- lgcp2$Lambda[,,1]
for(j in 2:(dim(lgcp2$Lambda)[3])){N <- N+lgcp2$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp2$xyt, cex=0.8, runtime=10, add=TRUE, prevalent="orange")




