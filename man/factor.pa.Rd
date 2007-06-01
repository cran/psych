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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}% __ONLY ONE__ keyword per line
