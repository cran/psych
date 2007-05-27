\name{eigen.loadings}
\alias{eigen.loadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert eigen vectors and eigen values to the more normal (for psychologists) component loadings}
\description{ The default procedures for principal component returns values not immediately equivalent to the loadings from a factor analysis.  eigen.loadings translates them into the more typical metric of eigen vectors multiplied by the squareroot of the eigenvalues.   This lets us find pseudo factor loadings if we have used princomp or princ.
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
"eigen.loadings" <- 
function (x) { #convert eigen vectors to loadings by unnormalizing them
               #used if using princomp or princ, not needed for principal
    return(x$vectors \%*\% sqrt(diag(x$values)))
    }
}
    
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
