
plot.stpp <- function(x, s.region=NULL, t.region=NULL, ...)
{
if (inherits(x,"stpp")==TRUE) 
	{ 
	  par(mfrow=c(1,2),pty="s")
	  if (is.null(s.region))	
	  plot(x[,1:2],pch=20,main="xy-locations")
	  else
		{
		  polymap(s.region,xlab="x",ylab="y")
		  points(x[,1:2],pch=20)
		  title("xy-locations")	 
		}
	  plot(x[,3],cumsum(x[,3]),type="l",xlab="t",ylab="",main="cumulative number",las=1,xlim=t.region)
 	}
}
getS3method("plot", "stpp", optional = FALSE)

