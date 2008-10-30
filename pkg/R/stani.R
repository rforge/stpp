
.stan3d.redraw <- function(o) {
  ## switch off redraws
  par3d(skipRedraw=TRUE)
  np=dim(o$xyt)[1]

  ## compute new states
  tin = rep(1,np)
  tin[o$xyt[,3]>(o$t-o$width)]=2
  tin[o$xyt[,3]>o$t]=3

  ## which points have changed state since last time?
  changed = tin != o$xyt[,5]

  ## remove any that changed
  if(any(changed)){
    rgl.pop(id=o$xyt[changed,4])
  }

  ## now add them back in their correct state:
  for(i in (1:np)[changed]){
    material3d(o$states[[tin[i]]])
    o$xyt[i,4]=spheres3d(x=o$xyt[i,1],y=o$xyt[i,2],z=o$xyt[i,3],radius=o$states[[tin[i]]]$radius)
    o$xyt[i,5]=tin[i]
  }
  ## start drawing again:
  par3d(skipRedraw=FALSE)

  ## update the slider positions
  .store(o)
  
  ## and return the modified object:
  return(o)
}

.rp.stan3d <- function(xyt,tlim,twid,states) {
  t=tlim[1];width=twid
  e=new.env()
  stan.panel  <- rp.control(title="space-time animation",
                            xyt=xyt, t=tlim[1], width=twid,
                            states=states,
                            e=e
                            )
  rp.slider(stan.panel, t, title = "time", from=tlim[1], to=tlim[2], action = .stan3d.redraw,showvalue=TRUE)
  rp.slider(stan.panel, width, title = "window", from=0, to=diff(tlim), action = .stan3d.redraw,showvalue=TRUE)
  rp.button(stan.panel,action=function(p){par3d(userMatrix = rotationMatrix(0, 1,0,0));return(p)},title="align time axis")
  rp.button(stan.panel,action=.store,title="quit",quitbutton=TRUE)
  rp.do(stan.panel, .stan3d.redraw)
  rp.block(stan.panel)
  return(e)
}

.store = function(panel){
 assign("t",panel$t,env=panel$e)
 assign("width",panel$width,env=panel$e)
 return(panel)
}

stani <- function(xyt,tlim=range(xyt[,3],na.rm=TRUE),twid=diff(tlim)/20,persist=FALSE,states,bgpoly,bgframe=TRUE,bgimage,bgcol=gray(seq(0,1,len=12))){
  require(rgl)
  require(rpanel)
  if(missing(states)){
    ## default colouring scheme:
    states=list(
      s1=list(col="blue",radius=1/80,alpha=0.5,lit=FALSE),
      s2=list(col="red",radius=1/30,alpha=0.5,lit=FALSE),
      ## still-to-come points are invisible (alpha=0)
      s3=list(col="yellow",alpha=0.0,radius=1/80,lit=FALSE)
      )
    if(persist){
      states$s1=states$s2
    }
  }

  maxRadius = max(states$s1$radius,states$s2$radius,states$s3$radius)
  xr = range(xyt[,1],na.rm=TRUE)
  yr = range(xyt[,2],na.rm=TRUE)
  tr = range(xyt[,3],na.rm=TRUE)
  diag=sqrt(diff(xr)^2+diff(yr)^2+diff(tr)^2)

  states$s1$radius=states$s1$radius*diag
  states$s2$radius=states$s2$radius*diag
  states$s3$radius=states$s3$radius*diag
    
  .setPlot(xr[1],xr[2],yr[1],yr[2],tr[1],tr[2],maxRadius)

  if(!missing(bgpoly)){
    poly=rbind(bgpoly,bgpoly[1,])
    poly=cbind(poly,min(tr))
    lines3d(poly,size=2.0)
    if(bgframe){
      ci=chull(bgpoly)
      nci=length(ci)
      cpoints=bgpoly[ci,]
      cpoints2 = cpoints[rep(1:nci,rep(2,nci)),]
      cpoints2 = cbind(cpoints2,c(min(tr),max(tr)))
      segments3d(cpoints2)
      poly=cbind(poly[,1:2],max(tr))
      lines3d(poly,size=2.0)
    }
  }

  if(!missing(bgimage)){
    .setBG(bgimage,min(tr),col=bgcol)
  }
  
  xyt=data.frame(xyt)
  xyt$id=NA
  ## initially all points will need redrawing:
  xyt$state=-1

  
  ## these points will get redrawn immediately... probably a better way to do this:
  for(i in 1:(dim(xyt)[1])){
    xyt[i,4]=points3d(xyt[,1],xyt[,2],xyt[,3],alpha=0.0)
  }
  env = .rp.stan3d(xyt,tlim,twid,states)
  ret=list()
  for(n in ls(env)){
    ret[[n]]=get(n,env=env)
  }
  return(ret)
}

.ranger <- function(x,margin=0.2){
  lim=range(x,na.rm=TRUE)
  lim=lim+c(-margin,margin)*diff(lim)
  return(lim)
}

.setPlot=function(xmin,xmax,ymin,ymax,tmin,tmax,radius=1/20){
  require(rgl)
  diag=sqrt((xmax-xmin)^2+(ymax-ymin)^2+(tmax-tmin)^2)
  xr=c(xmax,xmin)+c(radius*(xmax-xmin),-radius*(xmax-xmin))*10
  yr=c(ymax,ymin)+c(radius*(ymax-ymin),-radius*(ymax-ymin))*10
  

  plot3d(xr,yr,c(tmin,tmax),type="n",col="red",box=FALSE,axes=TRUE,xlab="x",ylab="y",zlab="t")
  axis3d('x-')
  axis3d('y-')
  axis3d('z-')
  par3d(FOV=1)
  AR=(xmax-xmin)/(ymax-ymin)
  aspect3d(AR,1,1)
  par3d(userMatrix = rotationMatrix(0, 1,0,0))

}

.setBG=function(xyz,zplane,col=heat.colors(12)){
  cols = col[as.integer(cut(xyz$z,length(col)))]
  surface3d(xyz$x,xyz$y,rep(zplane,prod(dim(xyz$z))),col=cols,lit=FALSE)
}
