"multi.hist" <-
function(x,nrow=NULL,ncol=NULL,density=TRUE,freq=FALSE,bcol="white",dcol=c("black","black"),dlty=c("dashed","dotted"),main="Histogram, Density, and Normal Fit",...) {
if((!is.matrix(x)) & (!is.data.frame(x))) {nvar <- 1
    x <- as.matrix(x,ncol=1) } else {
    x <- as.data.frame(x) 
nvar <- dim(x)[2] } #number of variables
if (length(dcol)<2) dcol <- c(dcol,dcol) 
     if(!density & (main == "Histogram, Density, and Normal Fit")) main = "Histogram" 
     nsize=ceiling(sqrt(nvar))   #size of graphic
     if(is.null(nrow) ) {nrow <- nsize} else {ncol <- nvar/nrow}
     if(is.null(ncol)) {ncol  <- ceiling(nvar/nsize )} else {nrow <- nvar/ncol}
     
     old.par <- par(no.readonly = TRUE) # all par settings which can be changed
     par(mfrow=c(nrow,ncol))       #set new graphic parameters
     for (i in 1:nvar) {
    	 xlab=names(x)[i]                #get the names for the variables
    	if(density) {histo.density(x[,i],xlab=xlab,main=main,freq=freq,bcol,dcol=dcol,dlty=dlty,...)} else {
    				hist(x[,i],main=main,xlab=xlab,freq=freq,bcol,dcol=dcol,dlty=dlty,...)}
    }  #draw the histograms for each variable
     on.exit(par(old.par))   #set the graphic parameters back to the original
     }

"histo.density" <- 
function(x,main="Histogram, Density, and Normal Fit",freq=FALSE,xlab=NULL,bcol="white",dcol=c("black","black"),dlty=c("dashed","dotted"),...) {
h <-  hist(x,plot=FALSE)
m1 <- mean(x,na.rm=TRUE)
s1 <- sd(x,na.rm=TRUE)
d <- density(x,na.rm=TRUE)

 if(freq) {ymax <- max(h$count)} else {ymax <- max(h$density)}
 dmax <- max(d$y)
ymax <- max(ymax,dmax)
plot(h,freq=freq,ylim=c(0,ymax*1.2),main=main,xlab=xlab,col=bcol,...)
if(!freq) {lines(d,lty=dlty[1],col=dcol[1],...)
curve(dnorm(x,m1,s1),add=TRUE,lty=dlty[2],col=dcol[2],...)} else {
lines(d$x,lty=dlty[1],col=dcol[1],...)}
      
}