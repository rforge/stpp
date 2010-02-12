is.3dpoints <- function (p) 
{
    is <- FALSE
    	if (is.array(p)) 
      	  if (length(dim(p)) == 2) 
            	if (dim(p)[2] >= 3) 
                	is <- TRUE
    is
}
