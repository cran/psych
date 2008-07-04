"multi.hist" <-
function(x,nrow=NULL,ncol=NULL,density=TRUE,main="Histogram, Density, and Normal Fit") {
nvar <- dim(x)[2]  #number of variables
     if(!density & (main == "Histogram, Density, and Normal Fit")) main = "Histogram" 
     nsize=ceiling(sqrt(nvar))   #size of graphic
     if(is.null(nrow) ) nrow <-nsize
     if(is.null(ncol)) ncol  <- ceiling(nvar/nsize )
     
     old.par <- par(no.readonly = TRUE) # all par settings which can be changed
     par(mfrow=c(nrow,ncol))       #set new graphic parameters
     for (i in 1:nvar) {
    	 xlab=names(x)[i]                #get the names for the variables
    	if(density) {histo.density(x[,i],xlab=xlab,main=main)} else {
    				hist(x[,i],main=main,xlab=xlab)}
    }  #draw the histograms for each variable
     on.exit(par(old.par))   #set the graphic parameters back to the original
     }

"histo.density" <- 
function(x,main="Histogram, Density, and Normal Fit",xlab=NULL) {
h <-  hist(x,plot=FALSE)
ymax <- max(h$density)
m1 <- mean(x,na.rm=TRUE)
s1 <- sd(x,na.rm=TRUE)
d <- density(x,na.rm=TRUE)
ymax <- max(h$density)
dmax <- max(d$y)
ymax <- max(ymax,dmax)
plot(h,freq=FALSE,ylim=c(0,ymax*1.2),main=main,xlab=xlab)
lines(d,lty="dashed")
curve(dnorm(x,m1,s1),add=TRUE,lty="dotted")
}