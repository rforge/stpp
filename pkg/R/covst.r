covst <- function(dist,times,separable=TRUE,model,param=c(1,1,1,1,1,2),sigma2=1,scale=c(1,1),plot=TRUE,nlevels=10)
{

  nt <- length(times)
  np <- length(dist)

  model <- set.cov(separable,model,param,var.grf)
  
  gs <- array(0, dim = c(np,nt))
  storage.mode(gs) <- "double"

#  dyn.load("/home/gabriel/functions/SimulSTPP/libF/covst.dll")

  gs <- .Fortran("covst",
                 (gs),
                 as.double(dist),
                 as.integer(np),
                 as.double(times),
                 as.integer(nt),
                 as.integer(model),
                 as.double(param),
                 as.double(sigma2),
                 as.double(scale))[[1]]

  if (plot==TRUE)
    {
#      image(dist,times,gs,col=grey(((10*max(length(times),length(dist))):1)/(10*max(length(times),length(dist)))),xlab="h",ylab="t",cex.axis=1.5,cex.lab=2,font=2)
      image(dist,times,gs,col=grey((1000:1)/1000),xlab="h",ylab="t",cex.axis=1.5,cex.lab=2,font=2)
      contour(dist,times,gs,add=T,col=4,labcex=1.5,nlevels=nlevels)
    }
  
  return(gs)

}



