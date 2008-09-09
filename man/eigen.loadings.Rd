\name{eigen.loadings}
\alias{eigen.loadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert eigen vectors and eigen values to the more normal (for psychologists) component loadings}
\description{ The default procedures for principal component returns values not immediately equivalent to the loadings from a factor analysis.  eigen.loadings translates them into the more typical metric of eigen vectors multiplied by the squareroot of the eigenvalues.   This lets us find pseudo factor loadings if we have used princomp  or eigen. \cr
If we use \code{\link{principal}} to do our principal components analysis, then we do not need this routine.}
\usage{
eigen.loadings(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the output from eigen or a list of class princomp derived from princomp}
}

\value{
  A matrix of Principal Component loadings more typical for what is expected in psychometrics.  That is, they are scaled by the square root of the eigenvalues.
}

\author{
 \email{  revelle@northwestern.edu           } \cr
 \url{http://personality-project.org/revelle.html}}
\note{Useful for SAPA analyses}
\examples{
x <- eigen(Harman74.cor$cov)
x$vectors[1:8,1:4]  #as they appear from eigen
y <- princomp(covmat=Harman74.cor$cov) 
y$loadings[1:8,1:4] #as they appear from princomp
eigen.loadings(x)[1:8,1:4] # rescaled by the eigen values
}
    
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
