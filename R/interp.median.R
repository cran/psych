"interp.median"  <- 
function(x,w=1,na.rm=TRUE) {
  im <- interp.quantiles(x,q=.5,w,na.rm=na.rm)
  return(im)} 
  


"interp.quantiles"  <- 
function(x,q=.5,w=1,na.rm=TRUE) {
  if (!(q>0) | !(q<1) ) {stop("quantiles most be greater than 0 and less than 1 q = ",q)}
   if(is.vector(x)) {im <- interp.q(x,q,w,na.rm=na.rm) } else {
 if((is.matrix(x) | is.data.frame(x)) ){
     n <- dim(x)[2]
     im <- matrix(NA,ncol=n)
     for (i in 1:n) {im[i] <- interp.q(x[,i],q,w=w,na.rm=na.rm)}
     colnames(im) <- colnames(x)
     } else {stop('The data must be either a vector, a matrix, or a data.frame')}    
       return(im)
     }}
 
 "interp.q" <- 
    function(x,q=.5,w=1,na.rm=TRUE) {
   
    if(na.rm) { x <- x[!is.na(x)]}
      n <- length(x)
   n2 <- (n+1)*q
   o <- order(x)
   x <-   x[o] 
   ml <- x[floor(n2)]
   mh <- x[ceiling(n2)]
   m <- (mh+ml)/2
   xb <- sum(x<m)
   xa <- sum(x>m)
   am <- n - xa - xb
    if(am >1) {
       im <- m -.5 *w + w*(n*q - xb )/am    #even number at median           
       } else {im <- m }  #no ties    
       return(im)
    }
 
"interp.quart" <- 
   function(x,w=1,na.rm=TRUE) {
    q <- c(.25,.5,.75)
     if(na.rm) { x <- x[!is.na(x)]}
   n <- length(x)
   n2 <- (n+1)*q
   N<- n*q
   o <- order(x)
   x <-   x[o] 
   ml <- x[floor(n2)]
   mh <- x[ceiling(n2)]
   m <- (mh+ml)/2
   im<-  xa <- xb <- rep(NA,3)
   for (i in 1:3)  {xb[i] <- sum(x <m[i])
       xa[i] <- sum(x > m[i]) }
   am <- n - xa - xb
   for (i in 1:3) {if(am[i] >1) {
            im[i] <- m[i] - .5*w + w*(N[i]-xb[i])/am[i]} else im[i] <- m[i]}
     return(im)} 
 
   
"interp.quartiles"  <- 
function(x,w=1,na.rm=TRUE) {
   q <- c(.25,.5,.75)
    if(is.vector(x)) {im <- interp.quart(x,w,na.rm=na.rm) 
         names(im) <- c("Q1","Median","Q3") } else {
    nvar <- dim(x)[2]
    im <- matrix(NA,ncol=3,nrow=nvar)
    for (i in 1:nvar )  {
    im[i,] <- interp.quart(x[,i],w,na.rm=na.rm)}
   rownames(im) <- colnames(x)
   colnames(im) <- c("Q1","Median","Q3")}
  return(im)} 
 
  

"interp.values" <- 
function(x,w=1,na.rm=TRUE) {
n <- length(x)
tabx <- table(x)
cv <- as.numeric(names(tabx))
k <- 1
v <- x[order(x)]
for (i in 1:length(tabx)) {
 for (j in 1:tabx[i]) {
  v[k] <- i - .5* w + j/(tabx[i]+1) 
  k <- k+1 }
  }
return(v) }



"interp.boxplot" <- 
   function(x,w=1,na.rm=TRUE) {
   stats <- interp.quartiles(x,w,na.rm=na.rm)
 return(stats) }
 
 
 "interp.qplot.by" <- 
function(y,x,w=1,na.rm=TRUE,xlab="group",ylab="dependent",ylim=NULL,arrow.len=.05,typ="b",add=FALSE,...) { 
  z <- by(y,x,interp.quartiles)
  zname <- names(z)
  z <- matrix(unlist(z),ncol=3,byrow=TRUE)
  rownames(z) <- zname
  colnames(z) <- c("Q1","Median","Q3")
  xv <- as.numeric(zname)
  ymin <- min(y)
  ymax <- max(y)
  if(is.null(ylim))  {ylim <- c(ymin,ymax)}
  if(add) {
  points(xv,z[,2],typ=typ,...)} else {
  plot(xv,z[,2],ylim = ylim,xlab=xlab,ylab=ylab,typ=typ,...)}
  lenx <- length(xv)
    for (i in 1:lenx)  {
          xlow <- z[i,1]
         xcen <- z[i,2]
    	 xup <- z[i,3]
    	
        arrows(xv[i],xlow,xv[i],xup,length=arrow.len, angle = 90, code=3, lty = NULL, lwd = par("lwd"), xpd = NULL,...) 
    	 }
  }
 



