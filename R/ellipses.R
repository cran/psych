"ellipses" <- 
function(x,y=NULL,add=FALSE,smooth=TRUE, lm=FALSE,data=TRUE,n=2,span=2/3, iter=3,col="red", xlab =NULL,ylab= NULL,...) {
#based upon John Fox's ellipse functions
  done=FALSE   #this is a kludge!
  segments=51
 
  if((is.matrix(x)) | (is.data.frame(x))) { 
       if (dim(x)[2] >2 ) { pairs.panels(x) 
                           done=TRUE
                           } else {
                               if(is.null(xlab)) xlab=colnames(x)[1]
                               if(is.null(ylab)) ylab=colnames(x)[2]
         						y <- x[,2]
        						x <- x[,1]
        						
                          } #dim ==2
         } else {
		      if((!is.vector(x)) | (!is.vector(y))) {stop("x and y must be vectors") } }
 if(!done){ xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
   r = (cor(x, y,use="pairwise"))
    if(is.null(xlab)) xlab=colnames(x)
   if(is.null(ylab)) ylab=colnames(y)
  angles <- (0:segments) * 2 * pi/segments
  
    unit.circle <- cbind(cos(angles), sin(angles))

  if (abs(r)>0 )theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
  
  shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
   ellipse <- unit.circle %*% shape 

   ellipse[,1] <- ellipse[,1]*xs + xm
   ellipse[,2] <- ellipse[,2]*ys + ym
   
if (add) {
        lines(ellipse,col=col, ...)
       if(data) { points(xm,ym,pch=20,cex=1.5,col=col)}
    
} else {plot(x,y,xlab=xlab,ylab=ylab,...)
     points(xm,ym,pch=20,cex=1.5,col=col)
    lines(ellipse, type = "l",col=col,...)}
   
   if(smooth) { ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col, ...)}
   if(lm) {
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col,lty="dashed", ...)}
    
 if(n>1) { ellipse <- unit.circle %*% shape  #draw another one
  
   ellipse[,1] <- ellipse[,1]*2*xs + xm
   ellipse[,2] <- ellipse[,2]*2*ys + ym
   lines(ellipse,col=col, ...)
     }}
  
}


  #then, just a fun function to draw Minkowski "circles"
  
"minkowski"<- function(r=2,add=FALSE,main=NULL,xl=1,yl=1) {
 segments=51
 x <- cos((0:segments) * pi/(2* segments))  #this spaces the points out to make a cleaner drawing
 min.circle <- cbind(x*xl, yl*((1-x^r)/(x^r+(1-x^r)))^(1/r))
if(add) {points(min.circle,type="l") } else plot(min.circle,ylim=c(-1,1),xlim=c(-1,1),typ="l", xlab="",ylab="",main=main)
 points(-min.circle,typ="l")
 points(-min.circle[,1],min.circle[,2],typ="l")
 points(min.circle[,1],-min.circle[,2],typ="l")
}
