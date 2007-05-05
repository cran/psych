\name{eigen.loadings}
\alias{eigen.loadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract eigen vectors, eigen values, show loadings}
\description{ Some SAPA matrices are ill formed.  This lets us find pseudo factor loadings if we have used princomp or princ.
If we use \code{\link{principal}} to do our principal components analysis, then we do not need this routine.}
\usage{
eigen.loadings(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix of loadings from princ or princomp}
}

\value{
  A matrix of PA loadings more typical for what is expected in psychometrics.
}

\author{
 \email{  revelle@northwestern.edu           } \cr
 \url{http://personality-project.org/revelle.html}}
\note{Useful for SAPA analyses}
\examples{

## The function is currently defined as
function (x) { #convert eigen vectors to loadings by unnormalizing them
               #used if using princomp or princ, not needed for principal
    n <- length(x$values)
    x$values[ x$values<0] <- 0
    fix<-sqrt(x$values)
    result<- x$vectors * rep(fix, each = n)}
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
