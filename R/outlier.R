#created 27/7/14
"outlier" <- 
 function(x,plot=TRUE,bad=5,na.rm=TRUE,xlab,ylab,...) {
 if(missing(xlab)) xlab <- expression("Quantiles of " * ~chi ^2)
 if(missing(ylab)) ylab <- expression("Mahalanobis " * D^2)
 rn <- rownames(x)
 nvar <- ncol(x)
 n.obs <- nrow(x)
 if(!is.matrix(x)) x <- as.matrix(x)
  nvar <- ncol(x)
 Sx <- cov(x,use="pairwise")
 Sx.inv <- solve(Sx)
# Mx <- colMeans(x,na.rm=na.rm)
# x <- sweep(x,2,Mx)
#x <- t(scale(t(x),scale=FALSE))
x <- scale(x,scale=FALSE)
D2 <- t(apply(x,1,function(xx) colSums(xx * Sx.inv,na.rm=TRUE)))
D2 <- rowSums(D2*x,na.rm=TRUE)
 names(D2) <- rn
 
 if(plot) {
  Chi2 <- qchisq(ppoints(n.obs), df =  nvar)
  qqplot(Chi2, D2,
       main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
                         " vs. quantiles of" * ~ chi[nvar]^2),xlab=xlab,ylab=ylab,...)
abline(0, 1, col = 'gray')
worst <- order(D2,decreasing=TRUE)
text(Chi2[n.obs:(n.obs-bad+1)],D2[worst[1:bad]],names(D2)[worst[1:bad]],pos=3,...)
}
 return(D2)
 }
 
