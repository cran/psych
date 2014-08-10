"violinBy" <- "densityBy" <- function(x,grp=NULL,grp.name=NULL,ylab="Observed",xlab="",main="Density plot",density=20,restrict=TRUE,xlim=NULL,add=FALSE,col=NULL,pch=20, ...) {
SCALE=.3  #how wide are the plots?

if(missing(col)) {col <- c("blue","red")}
nvar <- nvarg <- ncol(x)
if(!is.null(grp)) {
 if(!is.data.frame(grp) && !is.list(grp) && (length(grp) < NROW(x))) grp <- x[,grp]	
 Qnt <-  apply(x,2,function(xx) by(xx,grp,quantile,prob=c(0,1,.5,.25,.75),na.rm=TRUE))
meanX <- apply(x,2,function(xx) by(xx,grp,mean,na.rm=TRUE))
meanX <- matrix(unlist(meanX))
   Qnt <- matrix(unlist(Qnt),nrow=5)
   ngrp <- ncol(Qnt)/nvar
   nvarg <- ncol(Qnt)
} else {Qnt <- apply(x,2,quantile,prob=c(0,1,.5,.25,.75),na.rm=TRUE)
meanX <- apply(x,2,mean,na.rm=TRUE)}
minx <- Qnt[1,]
maxx <- Qnt[2,]
medx <- Qnt[3,]
Q25 <-  Qnt[4,]
Q75 <-  Qnt[5,]
rangex <- apply(x,2,range,na.rm=TRUE)
names <- colnames(x)

if(!is.null(grp)) {
if(missing(grp.name)) grp.name <- 1:ngrp
names <- paste(rep(names,each=ngrp),grp.name[1:ngrp],sep=" ")
     col <- rep(col,nvar* ngrp)}
d <- list(nvar)
if(missing(xlim)) xlim <- c(.5,nvarg+.5)

for (i in 1:nvar) {
if(!is.null(grp)) { if(restrict) {d[[i]] <- by(x[,i], grp ,function(xx) density(xx,na.rm=TRUE,from=rangex[1,i],to=rangex[2,i]))} else {
d[[i]] <- by(x[,i], grp ,function(xx) density(xx,na.rm=TRUE)) }} else {
if(restrict) {d[[i]] <- density(x[,i],na.rm=TRUE,from=minx[i],to=maxx[i])} else {
d[[i]] <- density(x[,i],na.rm=TRUE)} }
}

if(!add) {plot(meanX,ylim=c(min(minx),max(maxx)),xlim=xlim,axes=FALSE,xlab=xlab,ylab=ylab,main=main,pch=pch,...)
  axis(1,1:nvarg,names)
  axis(2)
  box()}
if(!is.null(grp)) d <- unlist(d,recursive=FALSE)
rev <- (length(d[[1]]$y):1) #this prevents a line down the middle
for(i in 1:nvarg) {
width <- SCALE/max(d[[i]]$y)
polygon(width*c(-d[[i]]$y,d[[i]]$y[rev])+i,c(d[[i]]$x,d[[i]]$x[rev]),density=density,col=col[i],...)
dmd <- max(which(d[[i]]$x<medx[i]))
d25 <- max(which(d[[i]]$x <= Q25[i]))
d75 <- max(which(d[[i]]$x <= Q75[i]))
segments(x0=width*d[[i]]$y[dmd] +i ,y0=d[[i]]$x[dmd],x1=-width*d[[i]]$y[dmd]+i,y1=d[[i]]$x[dmd],lwd=2)
segments(x0=width*d[[i]]$y[d25] +i ,y0=d[[i]]$x[d25],x1=-width*d[[i]]$y[d25]+i,y1=d[[i]]$x[d25])
segments(x0=width*d[[i]]$y[d75] +i ,y0=d[[i]]$x[d75],x1=-width*d[[i]]$y[d75]+i,y1=d[[i]]$x[d75])
}
 
 }
#created March 10, 2014 following a discussion of the advantage of showing distributional values
#modified March 22, 2014 to add the grouping variable
#basically just a violin plot. 


