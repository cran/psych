#revised Sept 16, 2013 to give control over the positon (pos) and size (cex) of the labels
#revised June 21, 2016 to allow naming of the points
 "cluster.plot" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Cluster plot",pch=18,pos,show.points=TRUE,choose=NULL,...) {
 if (!is.matrix(ic.results) ) {if (!is.null(class(ic.results)) )   {
  if(inherits(ic.results[1],"kmeans")) { load <- t(ic.results$centers) }  else {
      load <-ic.results$loadings} }} else {load <- ic.results}
      if(!is.null(choose)) load <- load[,choose,drop=FALSE]
nc <- dim(load)[2]
nvar <- dim(load)[1]

#defined locally, so as to be able to pass a parameter to it
  "panel.points" <- 
function (x, y,  pch = par("pch"), ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
   if(show.points) points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)
    text(x,y,vnames,...)
}


if(missing(pos)) pos <- rep(1,nvar)   #this allows more control over plotting
ch.col=c("black","blue","red","gray","black","blue","red","gray")
if (is.null(cluster)) {
cluster <- rep(nc+1,nvar)
cluster <- apply( abs(load)  ,1,which.max)
cluster[(apply(abs(load),1,max) < cut)] <- nc+1
}
if (nc > 2 ) {
vnames <- labels  #global variable
 pairs(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,lower.panel=panel.points,upper.panel=panel.points,...) }
 else {
if(show.points) { plot(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) } else {
  plot(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,type="n",...)
  pos=NULL}
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc < 3) text(load,labels,pos=pos,...)
 }
 
 


 
  "factor.plot" <-  
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title,jiggle=FALSE,amount=.02,pch=18,pos,show.points=TRUE,...) { #deprecated
fa.plot(ic.results,cluster=cluster,cut=cut,labels=labels,title=title,jiggle=jiggle,amount=amount,pch=pch,pos=pos,show.points=show.points,...)
}


 "fa.plot" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title,jiggle=FALSE,amount=.02,pch=18,pos,show.points=TRUE,choose=NULL,main=NULL,...) {
if(missing(title) ) { title="Plot"
              if (length(class(ic.results)) >1 )  {if (inherits(ic.results, "fa")) {title = "Factor Analysis"} else {
                     if (inherits(ic.results,"principal")) {title = "Principal Component Analysis"} 
                     } }
                     }
 if(missing(main)) {main<- title} else {title <- main}  #getting rid of confusion from years ago
 if (!is.matrix(ic.results)) {
        if (!is.null(class(ic.results))) {
       if(inherits(ic.results, "kmeans")) { load <- t(ic.results$centers) }  else {
      load <-ic.results$loadings} }} else {load <- ic.results}
      
      
   if(is.null(colnames(load))) colnames(load) <- paste("F",1:ncol(load),sep="")
if(!is.null(choose)) load <- load[,choose,drop=FALSE]
nc <- dim(load)[2]
nvar <- dim(load)[1]
if(missing(pos)) pos <- rep(1,nvar)   #this allows more control over plotting
ch.col=c("black","blue","red","gray","black","blue","red","gray")
if (is.null(cluster)) {
cluster <- rep(nc+1,nvar)
cluster <- apply( abs(load)  ,1,which.max)
cluster[(apply(abs(load),1,max) < cut)] <- nc+1
}
#define this function inside the bigger function so that vnames is globabl to it
 "panel.points" <- 
function (x, y,  pch = par("pch"), ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
   if(show.points) points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)

    text(x,y,vnames,...)
}
if(jiggle) load <- jitter(load,amount=amount)
if (nc > 2 ) {
vnames <- labels  #
 pairs(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],lower.panel=panel.points,upper.panel = panel.points,main=title,...) }
 else {
if(show.points) { plot(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) } else {
 
     plot(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,type="n",...)
     pos=NULL}
 abline(h=0)
 abline(v=0)
 if(is.null(labels)) labels <- paste(1:nvar)
  text(load,labels,pos=pos,...)}
 }
 
