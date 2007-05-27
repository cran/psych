\name{factor.pa}
\alias{factor.pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Principal Axis Factor Analysis}
\description{Among the many ways to do factor analysis, one of the most conventional is principal axes.  An eigen value decomposition of a correlation matrix is done and then the communalities for each variable are estimated by the first n factors. These communalities are entered onto the diagonal and the procedure is repeated until the sum(diag(r)) does not vary.  For well behaved matrices, maximum likelihood factor analysis (factanal) is probably preferred.  
}
\usage{
factor.pa(r, nfactors=1, residuals = FALSE, rotate = "varimax", min.err = 0.001, digits = 2, max.iter = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{A correlation matrix or a raw data matrix. If raw data, the correlation matrix will be found using pairwise deletion. }
  \item{nfactors}{ Number of factors to extract, default is 1 }
  \item{residuals}{Should the residual matrix be shown }
  \item{rotate}{"none", "varimax", "promax" }
  \item{min.err}{Iterate until the change in communalities is less than min.err}
  \item{digits}{ How many digits of output should be returned}
  \item{max.iter}{Maximum number of iterations for convergence }
}
\details{Factor analysis is an attempt to approximate a correlation or covariance matrix with one of lesser rank.  The basic model is that \eqn{nRn \approx nFk kFn'} where k is much less than n. There are many ways to do factor analysis, and maximum likelihood procedures are probably the most preferred (see factanal).  

Principal axes factor analysis has a long history in exploratory analysis and is a straightforward procedure.  Successive eigen value decompositions are done on a correlation matrix with the diagonal replaced with  diag (FF') until sum(diag(FF')) does not change (very much).  The current limit of max.iter =50 seems to work for most problems, but the Holzinger-Harmon 24 variable problem needs about 203 iterations to converge for a 5 factor solution.   

Principal axes may be used in cases when maximum likelihood solutions fail to converge.  
}
\value{
  If it is a LIST, use
  \item{values }{Eigen values of the final solution}
  \item{loadings }{An item by factor loading matrix.  Suitable for use in other programs (e.g., GPA rotation or factor2cluster.}
  \item{fit}{How well does the factor model reproduce the correlation matrix. (See VSS, ICLUST, and principal for this fit statistic.}
  \item{communality}{The history of the communality estimates.  Probably only useful for teaching what happens in the process of iterative fitting.}
}
\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.
}
\author{ William Revelle }

\seealso{ \code{\link{principal}}, \code{\link{VSS}}, \code{\link{ICLUST}}}
\examples{
#using the Harmon 24 mental tests, compare a principal factor with a principal components solution
pc <- principal(Harman74.cor$cov,4,rotate=TRUE)
pa <- factor.pa(Harman74.cor$cov,4,rotate="varimax")
round(factor.congruence(pc,pa),2)

#then compare with a maximum likelihood solution using factanal
mle <- factanal(x,4,covmat=Harman74.cor$cov)
round(factor.congruence(mle,pa),2)
#note that the order of factors and the sign of some of factors differ 

## The function is currently defined as
"factor.pa" <- 
function(r,nfactors=1,residuals=FALSE,rotate="varimax",min.err = .001,digits=2,max.iter=50) {
    n <- dim(r)[1]
    if (n!=dim(r)[2]) r <- cor(r,use="pairwise") # if given a rectangular matrix, the find the correlations first
    if (!residuals) { result <- list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)}

    orig <- diag(r)
    r.mat <- r
    comm <- sum(diag(r.mat))
    err <- comm
     i <- 1
    comm.list <- list()
    while(err > min.err)    #iteratively replace the diagonal with our revised communality estimate
      {
        eigens <- eigen(r.mat)
         loadings <- eigens$vectors[,1:nfactors] %*% diag(sqrt(eigens$values[1:nfactors]))
         model <- loadings \%*\% t(loadings)
         new <- diag(loadings \%*\% t(loadings))
         comm1 <- sum(new)
         diag(r.mat) <- new
         err <- abs(comm-comm1)
         comm <- comm1
         comm.list[[i]] <- comm1
         i <- i + 1
         if(i > max.iter) {warning("maximum iteration exceeded")
                err <-0 }
       }

       #make each vector signed so that the maximum loading is positive
    if (nfactors >1) {
    		maxabs <- apply(apply(loadings,2,abs),2,which.max)
   			sign.max <- vector(mode="numeric",length=nfactors)
    		for (i in 1: nfactors) {sign.max[i] <- sign(loadings[maxabs[i],i])}
    		loadings <- loadings \%*\% diag(sign.max)
   		
    	} else {
    		mini <- min(loadings)
   			maxi <- max(loadings)
    		if (abs(mini) > maxi) {loadings <- -loadings }
    		loadings <- as.matrix(loadings)
    	} #sign of largest loading is positive
    	colnames(loadings) <- paste("PA",1:nfactors,sep='')
    rownames(loadings) <- rownames(r)

    if(rotate != "none") {if (nfactors ==1) {warning("Must have at least two factors to rotate")} else {
    	if (rotate=="varimax") { 
   			loadings <- varimax(loadings)$loadings } else { 
     	if (rotate=="promax") {loadings <- promax(loadings)$loadings } 
     	}
     }}
    if(nfactors<1) nfactors <- n

    residual<- r - loadings \%*\% t(loadings)
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
    if (residuals) {result$residual <- round(residual,digits)}
   r2.off <- r
   diag(r2.off) <- 0
   r2.off <- sum(r2.off^2)
    diag(residual) <- 0
    rstar.off <- sum(residual^2)
    result$fit <- round(1-rstar2/r2,digits)
    result$fitoff <- round(1-rstar.off/r2.off,digits)
    result$values <- round(eigens$values,digits)
    result$loadings <- round(loadings,digits)
    result$communality <- round(unlist(comm.list),digits)
 
    return(result) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}% __ONLY ONE__ keyword per line
