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
  \item{f}{ A factor model matrix or a list of class loadings}
}
\details{The basic factor equation is \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2}. Residuals are just  R* = R - F'F. The residuals should be (but in practice probably rarely are) examined to understand the adequacy of the factor analysis.  When doing Factor analysis or Principal Components analysis, one usually continues to extract factors/components until the residuals do not differ from those expected from a random matrix.
}
\value{ rstar is the residual correlation matrix.
}

\author{
Maintainer: William Revelle <revelle@northwestern.edu>
}

\seealso{\code{\link{factor.pa}},  \code{\link{principal}}, \code{\link{VSS}}, \code{\link{ICLUST}}}
\examples{
fa2 <- factor.pa(Harman74.cor$cov,2,rotate=TRUE)
 fa2resid <- factor.residuals(Harman74.cor$cov,fa2)
 fa2resid[1:4,1:4] #residuals with two factors extracted
 fa4 <- factor.pa(Harman74.cor$cov,4,rotate=TRUE)
 fa4resid <- factor.residuals(Harman74.cor$cov,fa4)
 fa4resid[1:4,1:4] #residuals with 4 factors extracted

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
