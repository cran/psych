"scatter.hist" <-
function(x,y=NULL,smooth=TRUE,ab=FALSE, correl=TRUE,density=TRUE,ellipse=TRUE,digits=2,method,cex.cor=1,title="Scatter plot + histograms",xlab=NULL,ylab=NULL,...) {
old.par <- par(no.readonly = TRUE) # save default 
if(missing(xlab)) {
if(!is.null(colnames(x))) {xlab=colnames(x)[1]
                           ylab=colnames(x)[2]} else {xlab="V1"
                                                      ylab="V2"} }
                                                  
if (is.null(y)) {y <- x[,2]
                 x <- x[,1]} else {if(!is.null(dim(x))) {x <- x[,1,drop=TRUE]    
                                   if(!is.null(colnames(y))) ylab <- colnames(y)
                                   if(!is.null(dim(y))) {y <- y[,1,drop=TRUE]} } 
                                   }
xhist <- hist(x,breaks=11,plot=FALSE)
yhist  <- hist(y,breaks=11,plot=FALSE)

xrange <- range(x,na.rm=TRUE)
yrange <- range(y,na.rm=TRUE)

nf <- layout(matrix(c(2,4,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) #locations to plot
par(mar=c(5,4,1,1))    #first plot is in location 1 
plot(x,y,xlim=xrange,ylim=yrange,xlab=xlab,ylab=ylab,...)
if(ab) abline(lm(y~x))
if(smooth) {
 ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
       # lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
lines(stats::lowess(x[ok],y[ok]),col="red")}
if(ellipse) {ellipses(x,y,add=TRUE)}

par(mar=c(0,4,2,0)) 
mp <- barplot(xhist$density, axes=FALSE, space=0) 
#xhist <- hist(x,breaks=11,plot=TRUE,freq=FALSE,axes=FALSE,col="grey",main="",ylab="")
 tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
  if(class(tryd) != "try-error") {
 d$x <- (mp[length(mp)] - mp[1]+1) * (d$x - min(xhist$breaks))/(max(xhist$breaks)-min(xhist$breaks))
  # d$y <- d$y/max(d$y)
if(density)   lines(d)}
title(title)

par(mar=c(5,0,0,2)) 
mp <-   barplot(yhist$density, axes=FALSE, space=0, horiz=TRUE) 
 tryd <- try( d <- density(y,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
  if(class(tryd) != "try-error") {
  temp <- d$y
 d$y <- (mp[length(mp)] - mp[1]+1) * (d$x - min(yhist$breaks))/(max(yhist$breaks)-min(yhist$breaks))
  d$x <- temp
 if(density)    lines(d)
   }
 par(mar=c(3,1,1,1)) 
 if(correl) {
plot(1,1,type="n",axes=FALSE)
#plot(x,y)
med.x <- median(x,na.rm=TRUE)
med.y <- median(y,na.rm=TRUE)
if(missing(method)) method <- "pearson"
 r = (cor(x, y,use="pairwise",method=method))
 txt <- format(c(r, 0.123456789), digits=digits)[1]
 if(missing(cex.cor)) {cex <- 0.75/strwidth(txt)} else {cex <- cex.cor}
  text(1,1, txt,cex=cex)}
par(old.par)
}
#version of March 7, 2011
#revised Sept 7, 2013 to include method option in cor