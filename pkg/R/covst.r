set.cov <- function(separable,model,param,sigma2)
  {
    mods <- 0
    modt <- 0
    mod <- 0    

    models <- c("exponential","cauchy","stable","wave","gneiting","cesare","matern","none")

    for(i in 1:length(model))
      {
        M <- which(models==model[i])
        if (length(M)==0) stop("the model is not implemented")
      }

    for (i in 1:length(unique(model)))
      {
        if (((isTRUE(separable)) & ((model[i]==models[5]) | (model[i]==models[6]))) | ((!(isTRUE(separable))) & ((model[i]==models[1]) | (model[i]==models[2]) | (model[i]==models[3]) | (model[i]==models[4]) | (model[i]==models[7])))) stop("'stcov' does not match with 'model'")
      }

    if (isTRUE(separable))
      {
        if ((length(model)!=1) & (length(model)!=2))
          stop("for separable covariance functions, 'model' must be of length 1 or 2")
        if (length(model)==1)
          {
            if (model=="none")
              {
                mods <- 0
                modt <- 0
              }
            if (model=="exponential")
              {
                mods <- 1
                modt <- 1
              }
            if (model=="stable")
              {
                mods <- 2
                if ((param[1] >2) | (param[1]<0)) stop("Stable model parameter must lie in [0,2]")
                modt <- 2
                if ((param[2] >2) | (param[2]<0)) stop("Stable model parameter must lie in [0,2]")
              }
            if (model=="cauchy")
              {
                mods <- 3
                if (param[1]<=0) stop("Cauchy model parameter must be strictly positive")
                modt <- 3
                if (param[2]<=0) stop("Cauchy model parameter must be strictly positive")
              }
            if (model=="wave")
              {
                mods <- 4
                modt <- 4
                }
            if (model=="matern") 
			{
                  mods <- 7
			if (param[2]<=0 | param[1]<=0) stop("Matérn model parameters must be strictly positive")
                  modt <- 7
			if (param[3]<=0 | param[4]<=0) stop("Matérn model parameters must be strictly positive")
                  }
          }
            if (length(model)==2)
              {
                if (model[1]=="none")
                    mods <- 0
                if (model[2]=="none")
                  modt <- 0
                if (model[1]=="exponential")
                    mods <- 1
                if (model[2]=="exponential")
                  modt <- 1
                if (model[1]=="stable")
                    {
                      mods <- 2
                      if ((param[1] >2) | (param[1]<0)) stop("Stable model parameter must lie in [0,2]")
                    }
                if (model[2]=="stable")
                  {
                    modt <- 2
                    if ((param[2] >2) | (param[2]<0)) stop("Stable model parameter must lie in [0,2]")
                  }
                if (model[1]=="cauchy")
                  {
                    mods <- 3
                    if (param[1]<=0) stop("Cauchy model parmaeter must be strictly positive")
                  }
                if (model[2]=="cauchy")
                  {
                    modt <- 3
                    if (param[2]<=0) stop("Cauchy model parameter must be strictly positive")
                  }
                if (model[1]=="wave")
                    mods <- 4
                if (model[2]=="wave")
                  modt <- 4
                if (model[1]=="matern")
			{
                  mods <- 7
			if (param[2]<=0 | param[1]<=0) stop("Matérn model parameters must be strictly positive")
                  }
		    if (model[2]=="matern")
			{
                  modt <- 7
			if (param[3]<=0 | param[4]<=0) stop("Matérn model parameters must be strictly positive")
                  }
              }
      }
    if (!(isTRUE(separable)))
      {
        if (length(model)!=1)
          stop("for non-separable covariance functions, 'model' must be of length 1")
        if (model=="gneiting")
          {
            mod <- 5
            if (param[6]<2) stop("for Gneiting's covariance function, the sixth parameter must be greater than 2")
            if ((param[3]<=0) | (param[3]>2)) stop("for Gneiting's covariance function, the third parameter must lie in (0,2]")
            if ((param[4]<=0) | (param[4]>1)) stop("for Gneiting's covariance function, the fourth parameter must lie in (0,1]")
            if ((param[5]!=1) & (param[5]!=2) & (param[5]!=3)) stop("for Gneiting's covariance function, the fifth parameter must be 1, 2 or 3")
            if ((param[2]!=1) & (param[2]!=2) & (param[2]!=3)) stop("for Gneiting's covariance function, the second parameter must be 1, 2 or 3")
            if ((param[2]==1) & ((param[1]<0) | (param[1]>2))) stop("for Gneiting's covariance function, if the second parameter equals 1, the first parameter must lie in [0,2]") 
            if ((param[2]==2) & (param[1]<=0)) stop("for Gneiting's covariance function, if the second parameter equals 2, the first parameter must be strictly positive")            
          }
        if (model=="cesare")
          {
            mod <- 6
            if (((param[1]>2) | (param[1]<1)) | ((param[2]>2) | (param[2]<1))) stop("for De Cesare's model, the first and second parameters must lie in [1,2]")
            if (param[3]<3/2) stop("for De Cesare's model, the third parameter must be greater than 3/2")
          }
      }

    return(model=c(mods,modt,mod))
  }

matern = function (d, scale = 1, alpha = 1, nu = 0.5) 
{
    if (any(d < 0)) 
        stop("distance argument must be nonnegative")
    d <- d * alpha
    d[d == 0] <- 1e-10
    k <- 1/((2^(nu - 1)) * gamma(nu))
    res <- scale * k * (d^nu) * besselK(d, nu)
    return(res)
}


covst <- function(dist,times,separable=TRUE,model,param=c(1,1,1,1,1,2),sigma2=1,scale=c(1,1),plot=TRUE,nlevels=10)
{

  nt <- length(times)
  np <- length(dist)

  model <- set.cov(separable,model,param,sigma2)
  
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
 if (model[1]==7) mods=matern(dist,scale=scale[1],alpha=param[2],nu=param[1])
 if (model[2]==7) modt=matern(times,scale=scale[2],alpha=param[4],nu=param[3])
 if (model[1]==7 & model[2]==7) gs=mods%*%t(modt)
 if (model[1]==7 & model[2]!=7) gs=mods*gs
 if (model[2]==7 & model[1]!=7) gs=gs*modt


  if (plot==TRUE)
    {
#      image(dist,times,gs,col=grey(((10*max(length(times),length(dist))):1)/(10*max(length(times),length(dist)))),xlab="h",ylab="t",cex.axis=1.5,cex.lab=2,font=2)
      image(dist,times,gs,col=grey((1000:1)/1000),xlab="h",ylab="t",cex.axis=1.5,cex.lab=2,font=2)
      contour(dist,times,gs,add=T,col=4,labcex=1.5,nlevels=nlevels)
    }
  
  return(gs)

}



