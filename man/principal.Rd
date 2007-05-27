\name{principal}
\alias{principal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Principal components analysis }
\description{Does an eigen value decomposition and returns eigen values, loadings, and degree of fit for a specified number of components.  Basically just is doing a principal components for n principal components.  Can show the residual correlations as well. The quality of reduction in the squared correlations is reported by comparing residual correlations to original correlations. Unlike princomp, this returns a subset of just the best nfactors. The eigen vectors are rescaled by the sqrt of the eigen values to produce the component loadings more typical in factor analysis.
}
\usage{
principal(r, nfactors = 0, residuals = FALSE,rotate=FALSE, digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a correlation matrix }
  \item{nfactors}{Number of components to extract }
  \item{residuals}{ FALSE, do not show residuals, TRUE, report residuals }
  \item{rotate} {if TRUE, do varimax rotation, else report unrotated results}
  \item{digits }{digits =2    Accuracy of answers as well as display}
}
\details{Useful for those cases where the correlation matrix is improper (perhaps because of SAPA techniques).
}
\value{
  
  \item{values }{Eigen Values of all components -- useful for a scree plot}
  \item{loadings }{A standard loading matrix}
  \item{fit }{Fit of the model to the correlation matrix }
  \item{residual }{Residual matrix -- if requested}
}

\author{ William Revelle}

\seealso{\code{\link{VSS}},\code{\link{factor2cluster}},\code{\link{factor.pa}}, \code{\link{factor.congruence}}}
\examples{
#Four principal components of the Harmon 24 variable problem
#compare to a four factor principal axes solution using factor.congruence
pc <- principal(Harman74.cor$cov,4,rotate=TRUE)
pa <- factor.pa(Harman74.cor$cov,4,rotate=TRUE)
round(factor.congruence(pc,pa),2)

## The function is currently defined as
function(r,nfactors=0,residuals=FALSE,rotate=FALSE,digits=2) {
   n <- dim(r)[1]
  
   if (n!=dim(r)[2]) r <- cor(r,use="pairwise") # if given a rectangular matrix, the find the correlations first
    
    if (!residuals) { result <- list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)}
    eigens <- eigen(r)    #call the eigen value decomposition routine
    result$values <- round(eigens$values,digits)
   if(nfactors >0) {loadings <- eigen.loadings(eigens)[,1:nfactors] } else {loadings <- eigen.loadings(eigens)
   nfactors <- n}
   	
   	if(nfactors >1) {
    	maxabs <- apply(apply(loadings,2,abs),2,which.max)
    	sign.max <- vector(mode="numeric",length=nfactors)
    	for (i in 1: nfactors) {sign.max[i] <- sign(loadings[maxabs[i],i])}
  
  		loadings <- loadings %*% diag(sign.max)  #sign of largest loading is positive 
      	} else {if(nfactors >0) {mini <- min(loadings)
   			maxi <- max(loadings)
    		if (abs(mini) > maxi) {loadings <- -loadings }
    		loadings <- as.matrix(loadings)
    	 }}
    	
    colnames(loadings) <- paste("PC",1:nfactors,sep='')
    rownames(loadings) <- rownames(r)


    if(rotate) {loadings <- varimax(loadings)$loadings }  
   
    residual<- factor.residuals(r,loadings)
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
    if (residuals) {result$residual <-residual}
    result$fit <- round(1-rstar2/r2,digits)
    result$loadings <- round(loadings,digits)
    return(result)
   }
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
