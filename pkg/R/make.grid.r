make.grid <- function(nx,ny,poly)
  {
  #
  # Generate a rectangular grid of points in a polygon
  #
  # Requires Splancs package.
  #
  # Arguments:
  #
  #   nx, ny: grid dimensions.
  #
  #   poly: two columns matrix specifying polygonal region containing
  #         all data locations. If poly is missing, the unit square is
  #         considered.
  #
  # Value:
  #
  #    x, y: numeric vectors giving the coordinates of the points
  #          of the rectangular grid,
  #
  #    X, Y: matrix containing x and y coordinates respectively,
  #
  #     pts: index of the grid points belonging to the polygon,
  #
  #  xinc, yinc: mesh size of the grid,
  #
  #   mask: nx*ny matrix of logicals. TRUE if the grid point belongs
  #         to the polygon.
  #  
  ##
  ## E. GABRIEL, 28/12/2005
  ##

#    library(splancs)

    if (missing(poly)) poly <- matrix(c(0,0,1,0,1,1,0,1),4,2,T)
    
    if ((nx < 2) || (ny < 2)) stop("the grid must be at least of size 2x2")
    
    xrang <- range(poly[, 1], na.rm = TRUE)
    yrang <- range(poly[, 2], na.rm = TRUE)
    xmin <- xrang[1]
    xmax <- xrang[2]
    ymin <- yrang[1]
    ymax <- yrang[2]
    
    xinc <- (xmax-xmin)/nx
    yinc <- (ymax-ymin)/ny
    
    xc <- xmin-xinc/2
    yc <- ymin-yinc/2
    xgrid <- rep(0,nx)
    ygrid <- rep(0,ny)
    xgrid[1] <- xc + xinc
    ygrid[1] <- yc + yinc

    for (i in 2:nx)
      {
        xgrid[i] <- xgrid[i-1]+xinc
      }
    for (i in 2:ny)
      {
        ygrid[i] <- ygrid[i-1]+yinc
      }

    yy <- matrix(xgrid,nx,ny)
    xx <- t(yy)
    yy <- matrix(ygrid,nx,ny)

    X <- as.vector(xx)
    Y <- as.vector(yy)

    poly <- rbind(poly,poly[1,])
    pts <- inpip(pts=cbind(X,Y),poly)

    X[pts] <- TRUE
    X[X!=TRUE] <- FALSE
    mask <- matrix(X,ncol=ny,nrow=nx,byrow=T)

    invisible(return(list(x=xgrid,y=ygrid,X=xx,Y=yy,pts=pts,xinc=xinc,yinc=yinc,mask=matrix(as.logical(mask),nx,ny))))
  }
