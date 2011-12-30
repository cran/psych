#call the generic biplot function to draw factor/component loadings and scores
"biplot.psych" <-
function(x, labels=NULL,cex=c(.75,1),main="Biplot",hist.col="cyan",xlim=c(-3,3),ylim=c(-3,3),...) {
if(is.null(x$scores)) stop("Biplot requires factor/component scores:")
oldop <- par(no.readonly = TRUE)  
if(is.null(labels)) {if(nrow(x$scores) >100) {labels = rep(".",dim(x$scores)[1] )} else {labels = rep("o",dim(x$scores)[1] )}}
n.dims <- dim(x$loadings)[2]

if (n.dims == 2) {
 
  biplot(x$scores[,1:2],x$loadings[,1:2],xlabs=labels, cex=cex,main=main,...) } else {
	op1 <- par(mfrow=c(n.dims,n.dims), mar=c(2,3,3,2))
	for (i in 1:n.dims) {
   		for (j in 1:n.dims){ 
  		   if(i==j) {h <- hist(x$scores[,i],freq=FALSE, main=colnames(x$loadings)[i],xlab="",ylab="",col=hist.col)
  		        breaks <- h$breaks; nB <- length(breaks)
   		 tryd <- try( d <- density(x$scores[,i],na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
     	if(class(tryd) != "try-error") {
    	 lines(d)}
  } else {
        biplot(x$scores[,c(j,i)],x$loadings[,c(j,i)],xlabs=labels,xlab="",ylab="",cex=cex,xlim=xlim,ylim=ylim,...)}
                           }
                        }
                        par(oldop) }
}
	
