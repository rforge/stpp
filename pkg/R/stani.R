
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

  ## and return the modified object:
  return(o)
}

.rp.stan3d <- function(xyt,tlim,twid,states) {
  t=tlim[1];width=twid
  stan.panel  <- rp.control(title="space-time animation",
                            xyt=xyt, t=tlim[1], width=twid,
                            states=states
                            )
  rp.slider(stan.panel, t, title = "time", from=tlim[1], to=tlim[2], action = .stan3d.redraw,showvalue=TRUE)
  rp.slider(stan.panel, width, title = "window", from=0, to=diff(tlim), action = .stan3d.redraw,showvalue=TRUE)
  rp.button(stan.panel,action=function(p){par3d(userMatrix = rotationMatrix(0, 1,0,0));return(p)},title="align time axis")
  rp.do(stan.panel, .stan3d.redraw)
  
}

stani=function(xyt,tlim=range(xyt[,3],na.rm=TRUE),twid=diff(tlim)/20,persist=FALSE,states){
  require(rgl)
  require(rpanel)
  if(missing(states)){
    ## default colouring scheme:
    states=list(
      s1=list(col="blue",radius=1/80,alpha=1.0),
      s2=list(col="red",radius=1/20,alpha=1.0),
      ## still-to-come points are invisible (alpha=0)
      s3=list(col="yellow",alpha=0.0,radius=1/80)
      )
    if(persist){
      states$s1=states$s2
    }
  }
  xyt=data.frame(xyt)
  xyt$id=NA
  ## initially all points will need redrawing:
  xyt$state=-1

  xlim=.ranger(xyt[,1])
  ylim=.ranger(xyt[,2])
  tlim=.ranger(xyt[,3])

  plot3d(xlim,ylim,tlim,xlab="",ylab="",zlab="",axes=FALSE)
  par3d(FOV=1)
  ## aspect ratio...
  AR=diff(xlim)/diff(ylim)
  aspect3d(AR,1,1)
  par3d(userMatrix = rotationMatrix(0, 1,0,0))
  ## these points will get redrawn immediately... probably a better way to do this:
  for(i in 1:(dim(xyt)[1])){
    xyt[i,4]=points3d(xyt[,1],xyt[,2],xyt[,3])
  }
  .rp.stan3d(xyt,tlim,twid,states)
}

.ranger=function(x,margin=0.2){
  lim=range(x,na.rm=TRUE)
  lim=lim+c(-margin,margin)*diff(lim)
  return(lim)
}

#n=100
#data=cbind(runif(n),runif(n),rnorm(n))

#stan3d(data,persist=FALSE)
