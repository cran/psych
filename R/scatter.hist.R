"scatter.hist" <- "scatterHist" <- 
function(x,y=NULL,smooth=TRUE,ab=FALSE, correl=TRUE,data=NULL, density=TRUE,means=TRUE, ellipse=TRUE,digits=2,method="pearson",cex.cor=1,cex.point=1,title="Scatter plot + density",
   xlab=NULL,ylab=NULL,smoother=FALSE,nrpoints=0,xlab.hist=NULL,ylab.hist=NULL,grid=FALSE,xlim=NULL,ylim=NULL,x.breaks=11,y.breaks=11,
   x.space=0,y.space=0,freq=TRUE,x.axes=TRUE,y.axes=TRUE,size=c(1,2),col=c("blue","red","black"),legend=NULL,alpha=.5,pch=21,
   show.d=TRUE,
   x.arrow=NULL,y.arrow=NULL,d.arrow=FALSE,cex.arrow=1,...) {
old.par <- par(no.readonly = TRUE) # save default 
 sb <- grp <- NULL
 main <- title
 if(inherits(x, "formula")) {  ps <- fparse(x)
   formula <- TRUE
   if(is.null(data)) stop("You must specify the data if you are using formula input") 
     y  <- ps$y
     x <- ps$x
     if(length(x) > 1) {
        byGroup <- TRUE
        freq <- FALSE
        grp <- x[-1]
        grp.name <- grp
        grp <- data[,grp,drop=FALSE]
        n.grp <- length(table(grp))
        
        x <- x[1]
        xy <-data[,c(x,y,grp.name),drop=FALSE]
         sb <- statsBy(xy,group=grp.name,cors=TRUE)
     } else {grp <- NULL
        byGroup <- FALSE
        Md <- NULL
        }
       
       if(is.null(xlab)) xlab <- x
       if(is.null(ylab)) ylab <- y 
       x <- data[,x,drop=FALSE]
       y <- data[,y,drop=FALSE]
      
   }

col <- adjustcolor(col,alpha.f =alpha)
n.obs <- sum(!is.na(x))
if(missing(xlab)) {
if(!is.null(colnames(x))) {xlab=colnames(x)[1]
                           ylab=colnames(x)[2]} else {xlab="V1"
                                                      ylab="V2"} }

if (is.null(y)) {y <- x[,2]
                 x <- x[,1]} else {if(!is.null(dim(x))) {x <- x[,1,drop=TRUE]    
                                  # if(!is.null(colnames(y))) ylab <- colnames(y)
                                   if(!is.null(dim(y))) {y <- y[,1,drop=TRUE]} } 
                 }   
                 
if(missing(ylab)) { ylab <- colnames(y) }                                                               
if(is.null(grp)) grp <-1 
#if((length(pch) ==1 )&&(pch ==".")) pch <- 45     #this is kludge because we add 1 to pch later

if(NROW(grp) > 1) { byGroup <- TRUE
   dx <- by(x,grp,function(xx) density(xx, adjust=1,na.rm=TRUE))
   dy <- by(y,grp,function(xx) density(xx, adjust=1,na.rm=TRUE))     #what does adjust do?  I am copying this from my densityBy function
   x.means <- by(x, grp, function (xx) mean(xx, na.rm=TRUE))    #could clean this up to use statsBy output
   y.means <- by(y, grp, function(xx) mean(xx, na.rm=TRUE) )
   x.sd <-  by(x, grp, function (xx) sd(xx, na.rm=TRUE))
   y.sd <-  by(y, grp, function (xx) sd(xx, na.rm=TRUE))
    x.n <-by(x,grp,function(xx) sum(!is.na(xx)))
    y.n   <-by(y,grp,function(xx) sum(!is.na(xx)))
    x.d <-  (x.means[2]-x.means[1])/sqrt(((x.sd[1]^2*(x.n[1]-1))+ x.sd[2]^2* (x.n[2]-1))/(x.n[1] + x.n[2]))   #changed to n-1 7/29/21 
   y.d <-  (y.means[2]-y.means[1])/sqrt(((y.sd[1]^2*(y.n[1]-1)+ y.sd[2]^2* (y.n[2]-2)))/(y.n[1] + y.n[2])) 
   
   R.inv <- solve(sb$rwg)
   dist <- c(x.d,y.d)
    Md <- sqrt(t(dist) %*% R.inv %*% dist)

   if(show.d) {x.arrow=round(x.d,2)
               y.arrow=round(y.d,2)
               
               }
  grp <- unlist(grp)
  dx.max <- dy.max <- -9999
  for(i in 1:length(dx) ) {
  
  dx.max <- max(dx.max, dx[[i]]$y)
  dy.max <- max(dy.max, dy[[i]]$y)}
} else {byGroup <- FALSE
      x.means <- y.means <-x.d <- y.d<- stats<- NA   #give them a value so we have something to report in the results
      }

                                   
xrange <- range(x,na.rm=TRUE)
yrange <- range(y,na.rm=TRUE)
if(missing(xlim)) xlim <- xrange
if(missing(ylim)) ylim <- yrange
 x.breaks <- seq(xlim[1],xlim[2],(xlim[2] - xlim[1])/x.breaks)
 y.breaks <- seq(ylim[1],ylim[2],(ylim[2] - ylim[1])/y.breaks)                                   
xhist <- hist(x,breaks=x.breaks,plot=FALSE)
yhist  <- hist(y,breaks=y.breaks,plot=FALSE)


nf <- layout(matrix(c(2,4,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) #locations to plot the scatter plot
        #first plot is in location 1 
        
#nf <- layout(matrix(c(3,4,6,1, 2,5),2,3,byrow=TRUE),c(3,3,1),c(1,3))      #two x, 1 y   
# nf <- layout(matrix(c(5,6,9,3,4,8,1, 2,7),3,3,byrow=TRUE),c(3,3,1),c(1,3,3)) #two x, 2 y
#layout.show(nf) 
#the first figure is the plot
par(mar=c(5,4,1,1))    #bottom, left, top, right 


if(smoother) {smoothScatter(x,y,nrpoints=nrpoints,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)} else {
 if((length(pch) ==1 )&&(pch ==".")){plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,bg=col[grp],pch=pch)} else {
     plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,col=col[grp],bg=col[grp],pch=pch+grp,cex=cex.point)} } 
       

if(grid) grid()
if(ab) abline(lm(y~x))
if(smooth) {
 ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
       # lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
lines(stats::lowess(x[ok],y[ok]),col="red")}
if(ellipse) {
  if(byGroup) { grp.values <- table(grp)    #fix to work with multiple (>20 groups)
  grp.names <- dimnames(grp.values)
  ngrp <- length(grp.names$grp)
   for (i in 1:ngrp)  {x.grp <- x[grp==grp.names[[1]][[i]]]
                       y.grp <- y[grp==grp.names[[1]][[i]]]
  #lowx <- x[grp==grp.names[[1]][[1]]]
  #lowy <- y[grp==grp.names[[1]][[1]]]
  #highx <-  x[grp==grp.names[[1]][[2]]]
 # highy <- y[grp==grp.names[[1]][[2]]]
  ellipses(x.grp,y.grp,smooth=FALSE,add=TRUE,n=1,size=size,col=col[i],...)

 }
if(d.arrow) dia.arrow(c(x.means[1],y.means[1]), c(x.means[2],y.means[2]),labels=round(Md,2),both=TRUE,cex=cex.arrow)
} else {
   #do something (but how?)
  # xy <- cbind(x,y,grp)
#temp <- by(xy,grp,function(xx) colMeans(xx,na.rm=TRUE))
#for (i in length(temp) ){
# points(temp[[i]][1:2],pch=25,cex=4)
  ellipses(x,y,add=TRUE,size=size )}
  }

 if(!missing(legend) & byGroup) {        #show the legend
 
 location <- c("topleft","topright","top","left","right")
 
 grp.names <- paste(grp.name,names(table(grp)))
 n.grp <- length(grp.names)
 leg.text <- grp.names

 legend(location[legend],legend=leg.text,col=col[1:n.grp],fill=col[1:n.grp],pch=(pch+1:n.grp), lty=c(1:n.grp))
}

#this next figure is the density of the x axis
par(mar=c(.75,4,2,1))  #the left and right here should match the left and right from above  (the location for the x density)

#now, if we have a grouping variable, then plot a barplot for each grp value


if(byGroup) {

plot(dx[[1]],main="",xlim=xrange,axes=FALSE,ylim =c(0,dx.max))
title(main,...)

 for (i in 1:length(dx)) {
   if(freq) {scal <- dx[[i]]$n
      dx[[i]]$y <- dx[[i]]$y*scal} 
    
     polygon(dx[[i]] ,col=col[i])
     #x.mean <- mean(dx[[i]]$x)      #this is based upon the density, not the data
     x.mean <- x.means[i] 
     y.mean <-  mean(dx[[i]]$y[256])  #the median
      y.mean <- dx[[i]]$y[which.max(dx[[i]]$x > x.means[i])]  
     segments(x.mean,0,x.mean,y.mean)
      }
 
  if(!is.null(x.arrow)) {    
  dia.arrow(c(x.means[1],.2*max(dx[[1]]$y)), c(x.means[2],.2*max(dx[[1]]$y)),labels=x.arrow,both=TRUE,cex=cex.arrow)}
  
} else {
if(freq) { mp <- barplot(xhist$counts, axes=x.axes, space=x.space,xlab=xlab.hist)} else { mp <- barplot(xhist$density, axes=x.axes, space=x.space,xlab=xlab.hist)}
 #xhist <- hist(x,breaks=11,plot=TRUE,freq=FALSE,axes=FALSE,col="grey",main="",ylab="")
 tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
 
  if(!inherits(tryd ,"try-error")) {
 d$x <- (mp[length(mp)] - mp[1]+1) * (d$x - min(xhist$breaks))/(max(xhist$breaks)-min(xhist$breaks))
  if(freq) d$y <- d$y * max(xhist$counts/xhist$density,na.rm=TRUE)
  if(density)   lines(d)}
  


}
 
# title(title,cex=cex.title)   #doesn't actually do anything

#now show the y axis densities
par(mar=c(5,0.5,1,2))        #the location for the y density 

if(byGroup) {
 
 #if(freq) { mp <- barplot(xhist$counts, axes=x.axes, space=x.space,xlab=xlab.hist,plot=TRUE)} else { mp <- barplot(xhist$density, axes=x.axes, space=x.space,xlab=xlab.hist,plot=TRUE)}
 plot(dy[[1]],main="",axes=FALSE,ylim=yrange,xlim = c(0,dy.max),xlab="Density")
 
 

 for (i in 1:length(dy)) {
  
    
     temp <- dy[[i]]$y      #swap the x and y columns
     dy[[i]]$y <- dy[[i]]$x
     dy[[i]]$x <- temp
      if(freq) {scal <- dy[[i]]$n
      dy[[i]]$y <- dy[[i]]$y*scal} 
      
      
     polygon(dy[[i]] ,col=col[i])
     x.mean <- y.means[i]
     # x.mean <- mean(dy[[i]]$y)
    y.mean <- dy[[i]]$x[which.max(dy[[i]]$y > y.means[i])]   #this adjusts for the problem that means of density do not match real means

     segments(0,x.mean,y.mean,x.mean)
      }
      
       if(!is.null(y.arrow)) {    
  dia.arrow(c(.2*max(dx[[1]]$y),y.means[1]), c(.2*max(dx[[1]]$y), y.means[2]),labels=y.arrow,both=TRUE,cex=cex.arrow)}

  
  
} else {


if(freq) {mp <-   barplot(yhist$counts, axes=y.axes, space=y.space, horiz=TRUE,ylab=ylab.hist) } else {mp <-   barplot(yhist$density, axes=y.axes, space=y.space, horiz=TRUE,ylab=ylab.hist)}
 tryd <- try( d <- density(y,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
  if(!inherits(tryd,"try-error")) {
  temp <- d$y
 d$y <- (mp[length(mp)] - mp[1]+1) * (d$x - min(yhist$breaks))/(max(yhist$breaks)-min(yhist$breaks))
  d$x <- temp
  if(freq) d$x <- d$x * max(yhist$counts/yhist$density,na.rm=TRUE)
 if(density)    lines(d)
   }
}



#this is the upper right hand square  -- show the correlation/Md
 par(mar=c(1,1,1,1))    #the  location for the correlations if we are not doing groups


 if(correl) {
plot(1,1,type="n",axes=FALSE)
#plot(x,y)
med.x <- median(x,na.rm=TRUE)
med.y <- median(y,na.rm=TRUE)
if(missing(method)) method <- "pearson"
 r = (cor(x, y,use="pairwise",method=method))
# txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste0("r = ", round(r,digits),"\n")
 if(missing(cex.cor)) {cex <- 0.75/strwidth(txt)} else {cex <- cex.cor}
  text(1,1, txt,cex=cex)
  if(!is.null(sb))text(1,.8,paste0("r wg =",round(sb$rwg[1,2],2) ),cex=.8*cex)} else {
  plot(1,1,type="n",axes=FALSE)
  if(!is.null(Md) ) {
  
  txt <- paste0("D = ",sprintf("%.2f",round(Md,2)),"\n")
  if(missing(cex.cor)) {cex <- 0.75/strwidth(txt)} else {cex <- cex.cor}
  text(1,1,txt,cex=cex)
  #text(1,.8,paste0("r wg =",round(sb$rwg[1,2],2) ),cex=cex)
  }
  }
  
  result <- list(x.means=x.means,y.means=y.means,x.d=x.d, y.d = y.d,stats=sb)
par(old.par)
invisible(result)
}
#version of March 7, 2011
#revised Sept 7, 2013 to include method option in cor
#revised October 3, 2020 to allow groups and to allow formula input
#revised June 6, 2021  to allow ellipses by groups
