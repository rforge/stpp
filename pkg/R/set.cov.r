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
        if (((isTRUE(separable)) && ((model[i]==models[5]) || (model[i]==models[6]))) || ((!(isTRUE(separable))) && ((model[i]==models[1]) || (model[i]==models[2]) || (model[i]==models[3]) || (model[i]==models[4]) || (model[i]==models[7])))) stop("'stcov' does not match with 'model'")
      }

    if (isTRUE(separable))
      {
        if ((length(model)!=1) && (length(model)!=2))
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
                if ((param[1] >2) || (param[1]<0)) stop("parameter of stable model must lie in [0,2]")
                modt <- 2
                if ((param[2] >2) || (param[2]<0)) stop("parameter of stable model must lie in [0,2]")
              }
            if (model=="cauchy")
              {
                mods <- 3
                if (param[1]<=0) stop("parameter of cauchy model must be strictly greater than 0")
                modt <- 3
                if (param[2]<=0) stop("parameter of cauchy model must be strictly greater than 0")
              }
            if (model=="wave")
              {
                mods <- 4
                modt <- 4
                }
            if (model=="matern") stop("You must specify another model for the temporal covariance")
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
                      if ((param[1] >2) || (param[1]<0)) stop("parameter of stable model must lie in [0,2]")
                    }
                if (model[2]=="stable")
                  {
                    modt <- 2
                    if ((param[2] >2) || (param[2]<0)) stop("parameter of stable model must lie in [0,2]")
                  }
                if (model[1]=="cauchy")
                  {
                    mods <- 3
                    if (param[1]<=0) stop("parameter of cauchy model must be strictly greater than 0")
                  }
                if (model[2]=="cauchy")
                  {
                    modt <- 3
                    if (param[2]<=0) stop("parameter of cauchy model must be strictly greater than 0")
                  }
                if (model[1]=="wave")
                    mods <- 4
                if (model[2]=="wave")
                  modt <- 4
                if (model[1]=="matern")
                  mods <- 7
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
            if ((param[3]<=0) || (param[3]>2)) stop("for Gneiting's covariance function, the third parameter must lie in (0,2]")
            if ((param[4]<=0) || (param[4]>1)) stop("for Gneiting's covariance function, the fourth parameter must lie in (0,1]")
            if ((param[5]!=1) && (param[5]!=2) && (param[5]!=3)) stop("for Gneiting's covariance function, the fifth parameter must be 1, 2 or 3")
            if ((param[2]!=1) && (param[2]!=2)  && (param[2]!=3)) stop("for Gneiting's covariance function, the second parameter must be 1, 2 or 3")
            if ((param[2]==1) && ((param[1]<0) || (param[1]>2))) stop("for Gneiting's covariance function, if the second parameter equals 1, the first parameter must lie in [0,2]") 
            if ((param[2]==2) && (param[1]<=0)) stop("for Gneiting's covariance function, if the second parameter equals 2, the first parameter must be strictly positive")            
          }
        if (model=="cesare")
          {
            mod <- 6
            if (((param[1]>2) || (param[1]<1)) || ((param[2]>2) || (param[2]<1))) stop("for De Cesare's model, the first and second parameters must lie in [1,2]")
            if (param[3]<3/2) stop("for De Cesare's model, the third parameter must be greater than 3/2")
          }
      }

    return(model=c(mods,modt,mod))
  }
