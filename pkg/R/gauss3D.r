gauss3D <- function(nx=100,ny=100,nt=100,xlim,ylim,tlim,separable=TRUE,model="exponential",param=c(1,1,1,1,1,2),scale=c(1,1),var.grf=1,mean.grf=0,exact=TRUE)
{

  N <- c(nx,ny,nt)

  mod <- set.cov(separable,model,param,var.grf)
  
  g <- floor(log(2*(N-1))/log(2))+1
  M <- 2^g

  count <- 1
  changG <- TRUE
  while(changG==TRUE)
    {
      L <- rep(-9999,M[1]*M[2]*M[3])  
      
      storage.mode(L) <- "double"
      
      res <- .Fortran("circ",
                      (L),
                      as.integer(M),
                      as.integer(N),
                      as.double(xlim),
                      as.double(ylim),
                      as.double(tlim),
                      as.integer(mod), 
                      as.double(param),
                      as.double(var.grf),
                      as.double(scale))[[1]]
 
      L <- array(res,dim=M)
      FTL <- fft(L)

      if (isTRUE(exact))
        {      
          if (min(Re(FTL))<0)
            {
              g <- g+1
              M <- 2^g
              changG <- TRUE
              count <- count+1
            }
          else
            changG <- FALSE
        }
      else
        {
          FTL[Re(FTL)<0] <- 0
          changG <- FALSE
        }

    }
  print(count)
  
  X <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
  X <- fft(X,inverse=TRUE)
  A <- sqrt(FTL)*X
  G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
  G <- G+mean.grf

## Remarks
# --------
#
#  The right way is:
#    X1 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#    X2 <- array(rnorm(M[1]*M[2]*M[3],0,1),dim=M)
#    X <- fft(X1+1i*X2,inverse=TRUE)
#    A <- sqrt(FTL)*X
#    G <- (Re(fft(A,inverse=FALSE))/(M[1]*M[2]*M[3]))[1:N[1],1:N[2],1:N[3]]
#  but it doesn't change the results and it is slightly faster. 
#  
# Taking the first elements of FTL and then computing G (changing M by N)
# is faster than taking the first elements of G, but it does not provide
# the same results
  
  return(G)
}

