\name{factor.residuals}
\alias{factor.residuals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ R* =  R- F F' }
\description{The basic factor or principal components model is that a correlation or covariance matrix may be reproduced by the product of a factor loading matrix times its transpose.  Find the residuals of the original minus the  reproduced matrix.  Used by \code{\link{factor.fit}}, \code{\link{VSS}}, \code{\link{ICLUST}}, etc.
}
\usage{
factor.residuals(r, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{ A correlation matrix }
  \item{f}{ A factor model matrix }
}
\details{The basic factor equation is \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2}. Residuals are just  R* = R - F'F. 
}
\value{ rstar is the residual correlation matrix.
}

\author{
Maintainer: William Revelle <revelle@northwestern.edu>
}

\seealso{\code{\link{factor.pa}},  \code{\link{principal}}, \code{\link{VSS}}, \code{\link{ICLUST}}}
\examples{

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
