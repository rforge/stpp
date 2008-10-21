animation <- function(xyt, s.region, t.region, runtime=1, incident="red", prevalent="pink3", fade=0.01, pch=19, cex=0.25, plot.s.region=T, scales=T, border.frac=0.05, add=F)
{ 
  #
  # Description:
  #   Animation of spatio-temporal point process data.
  #   Requires splancs library
  #
  # Arguments:
  #            xyt: data - matrix containing the (x,y,t)-coordinates
  #        runtime: approximate running time of animation, in seconds
  #       s.region: two-column matrix specifying polygonal
  #                 region containing all data-locations xyt[,1:2]
  #       t.region: interval containing all data-times xyt[,3]
  #       incident: colour in which incident point xyt[i,1:2] is
  #                 plotted at time xyt[i,3]
  #      prevalent: colour to which prevalent point xyt[1,1:2] fades
  #                 at time xyt[i,3]+fade
  #           fade: approximate time-delay for colour-change from incident
  #                 to prevalent, in data-time units
  #            pch: plotting symbol usedfor each point 
  #            cex: magnification of plotting symbol relative to standard size
  #  plot.s.region: if true, plot s.region as polygon
  #         scales: if true, plot X and Y axes with scales
  #    border.frac: extent of border of plotting region surounding s.region,
  #                 as fraction of ranges of X and Y
  #
  #


  if (missing(s.region)) s.region <- sbox(xyt[,1:2],xfrac=0.01,yfrac=0.01)
  if (missing(t.region)) t.region <- range(xyt[,3],na.rm=T)
  
#  par(pty="s",mfrow=c(1,1))
  ott<-order(xyt[,3])
  sxyt<-xyt[ott,]
  rangex<-range(s.region[,1])
  rangey<-range(s.region[,2])
  xlim<-c(rangex[1]-border.frac*(rangex[2]-rangex[1]),rangex[2]+border.frac*(rangex[2]-rangex[1]))
  ylim<-c(rangey[1]-border.frac*(rangey[2]-rangey[1]),rangey[2]+border.frac*(rangey[2]-rangey[1]))
  xy<-as.matrix(sxyt[,1:2])
  tt<-sxyt[,3]
  npts<-length(tt)
  T0 <- max(t.region)

  if (add==F)
    {
      if (scales==F)
        plot(xy[,1],xy[,2],type="n",xlim=xlim,ylim=ylim,xaxt="n",yaxt="n",bty="n",xlab=" ",ylab=" ")
      if (scales==T)
        plot(sxyt[,1],sxyt[,2],type="n",xlim=xlim,ylim=ylim,bty="n",xlab="X",ylab="Y")
      if (plot.s.region==T)
        polymap(as.points(s.region),add=T,lwd=2)
    }
  nplotted<-0
  tt.now<-0
  tt.fade<-0
  while (nplotted<npts)
    {
      i<-nplotted+1
      tt.gap<-tt[i]-tt.now
      junk<-Sys.sleep((tt.gap/T0)*runtime)
      n.fade<-sum(tt[1:i]<tt.fade)
      if (sum(n.fade)>0)
        points(xy[1:n.fade,1],xy[1:n.fade,2],col=prevalent,pch=19,cex=cex)
      points(xy[i,1],xy[i,2],col=incident,pch=19,cex=cex)
      nplotted<-i
      tt.now<-tt[i]
      tt.fade<-tt.now-fade
    }
}

