"scatter.hist" <- "scatterHist" <- 
function(x,y=NULL,smooth=TRUE,ab=FALSE, correl=TRUE,data=NULL, density=TRUE,ellipse=TRUE,digits=2,method,cex.cor=1,title="Scatter plot + histograms",
   xlab=NULL,ylab=NULL,smoother=FALSE,nrpoints=0,xlab.hist=NULL,ylab.hist=NULL,grid=FALSE,xlim=NULL,ylim=NULL,x.breaks=11,y.breaks=11,
   x.space=0,y.space=0,freq=TRUE,x.axes=TRUE,y.axes=TRUE,size=c(1,2),col=c("blue","red","black"),legend=NULL,alpha=.5,pch=21,...) {
old.par <- par(no.readonly = TRUE) # save default 
 grp <- NULL
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
     } else {grp <- NULL
        byGroup <- FALSE
        }
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



if(NROW(grp) > 1) { byGroup <- TRUE
   dx <- by(x,grp,function(xx) density(xx, adjust=1,na.rm=TRUE))
   dy <- by(y,grp,function(xx) density(xx, adjust=1,na.rm=TRUE))     #what does adjust do?  I am copying this from my densityBy function
  grp <- unlist(grp)
  dx.max <- dy.max <- -9999
  for(i in 1:length(dx) ) {
  
  dx.max <- max(dx.max, dx[[i]]$y)
  dy.max <- max(dy.max, dy[[i]]$y)}
} else {byGroup <- FALSE}

                                   
xrange <- range(x,na.rm=TRUE)
yrange <- range(y,na.rm=TRUE)
if(missing(xlim)) xlim <- xrange
if(missing(ylim)) ylim <- yrange
 x.breaks <- seq(xlim[1],xlim[2],(xlim[2] - xlim[1])/x.breaks)
 y.breaks <- seq(ylim[1],ylim[2],(ylim[2] - ylim[1])/y.breaks)                                   
xhist <- hist(x,breaks=x.breaks,plot=FALSE)
yhist  <- hist(y,breaks=y.breaks,plot=FALSE)


nf <- layout(matrix(c(2,4,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE) #locations to plot the scatter plot
par(mar=c(5,4,1,1))    #first plot is in location 1 
if(smoother) {smoothScatter(x,y,nrpoints=nrpoints,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)} else {plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,bg=col[grp],pch=pch+grp)}

if(grid) grid()
if(ab) abline(lm(y~x))
if(smooth) {
 ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
       # lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
lines(stats::lowess(x[ok],y[ok]),col="red")}
if(ellipse) {
  #if(byGroup) {
   #browser()

   #do something (but how?)
  # xy <- cbind(x,y,grp)
#temp <- by(xy,grp,function(xx) colMeans(xx,na.rm=TRUE))
#for (i in length(temp) ){
# points(temp[[i]][1:2],pch=25,cex=4)
  ellipses(x,y,add=TRUE,size=size)}

 if(!missing(legend) & byGroup) {        #show the legend
 
 location <- c("topleft","topright","top","left","right")
 
 grp.names <- paste(grp.name,names(table(grp)))
 n.grp <- length(grp.names)
 leg.text <- grp.names

 legend(location[legend],legend=leg.text,col=col[1:n.grp],fill=col[1:n.grp],pch=(pch+1:n.grp), lty=c(1:n.grp))
}

par(mar=c(.75,4,2,1))  #the left and right here should match the left and right from above  (the location for the x density)

#now, if we have a grouping variable, then plot a barplot for each grp value


if(byGroup) {
 
plot(dx[[1]],main=main,xlim=xrange,axes=FALSE,ylim =c(0,dx.max))

 for (i in 1:length(dx)) {
   if(freq) {scal <- dx[[i]]$n
      dx[[i]]$y <- dx[[i]]$y*scal} 
   
     polygon(dx[[i]] ,col=col[i])
      }
   
   
  
} else {
if(freq) { mp <- barplot(xhist$counts, axes=x.axes, space=x.space,xlab=xlab.hist)} else { mp <- barplot(xhist$density, axes=x.axes, space=x.space,xlab=xlab.hist)}
 #xhist <- hist(x,breaks=11,plot=TRUE,freq=FALSE,axes=FALSE,col="grey",main="",ylab="")
 tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
 
  if(!inherits(tryd ,"try-error")) {
 d$x <- (mp[length(mp)] - mp[1]+1) * (d$x - min(xhist$breaks))/(max(xhist$breaks)-min(xhist$breaks))
  if(freq) d$y <- d$y * max(xhist$counts/xhist$density,na.rm=TRUE)
  if(density)   lines(d)}

}
 
 title(title)


par(mar=c(5,0.5,1,2))        #the locaction for the y density 


if(byGroup) {
 
 #if(freq) { mp <- barplot(xhist$counts, axes=x.axes, space=x.space,xlab=xlab.hist,plot=TRUE)} else { mp <- barplot(xhist$density, axes=x.axes, space=x.space,xlab=xlab.hist,plot=TRUE)}
temp <- dy[[1]]$y
dy[[1]]$y <- dy[[1]]$x
dy[[1]]$x <- temp
plot(dy[[1]],main="",axes=FALSE,ylim=yrange,xlim = c(0,dy.max),xlab="Density")
polygon(dy[[1]], col=col[1])

 for (i in 2:length(dy)) {
   if(freq) {scal <- dy[[i]]$n
      dy[[i]]$y <- dy[[i]]$y*scal} 
    
     temp <- dy[[i]]$y
     dy[[i]]$y <- dy[[i]]$x
     dy[[i]]$x <- temp
     polygon(dy[[i]] ,col=col[i])
      }
   
   
  
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

 par(mar=c(1,1,1,1))    #the  location for the correlations
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
#revised October 3, 2020 to allow groups and to allow formula input
