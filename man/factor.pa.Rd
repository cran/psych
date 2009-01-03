\name{factor.pa}
\alias{factor.pa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Principal Axis Factor Analysis}
\description{Among the many ways to do factor analysis, one of the most conventional is principal axes.  An eigen value decomposition of a correlation matrix is done and then the communalities for each variable are estimated by the first n factors. These communalities are entered onto the diagonal and the procedure is repeated until the sum(diag(r)) does not vary.  For well behaved matrices, maximum likelihood factor analysis (factanal) is probably preferred.  
}
\usage{
factor.pa(r, nfactors=1, residuals = FALSE, rotate = "varimax",n.obs = NA,
scores = FALSE,SMC=TRUE, missing=FALSE,impute="median",min.err = 0.001, digits = 2, max.iter = 50,symmetric=TRUE,warnings=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{A correlation matrix or a raw data matrix. If raw data, the correlation matrix will be found using pairwise deletion. }
  \item{nfactors}{ Number of factors to extract, default is 1 }
  \item{residuals}{Should the residual matrix be shown }
  \item{rotate}{"none", "varimax", "promax" or "oblimin" are possible rotations of the solution.}
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics.}
  \item{scores}{If TRUE, estimate factor scores}
  \item{SMC}{Use squared multiple correlations (SMC=TRUE) or use 1  as initial communality estimate. Try using 1 if imaginary eigen values are reported. }
  \item{missing}{if scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean}
  \item{impute}{"median" or "mean" values are used to replace missing values}
  \item{min.err}{Iterate until the change in communalities is less than min.err}
  \item{digits}{ How many digits of output should be returned}
  \item{max.iter}{Maximum number of iterations for convergence }
  \item{symmetric}{symmetric=TRUE forces symmetry by just looking at the lower off diagonal values}
  \item{warnings}{warnings=TRUE => warn if number of factors is too many }
}
\details{Factor analysis is an attempt to approximate a correlation or covariance matrix with one of lesser rank.  The basic model is that \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2} where k is much less than n. There are many ways to do factor analysis, and maximum likelihood procedures are probably the most preferred (see \code{\link{factanal}} ).  The existence of uniquenesses is what distinguishes factor analysis from principal components analysis (e.g., \code{\link{principal}}).  

Principal axes factor analysis has a long history in exploratory analysis and is a straightforward procedure.  Successive eigen value decompositions are done on a correlation matrix with the diagonal replaced with  diag (FF') until sum(diag(FF')) does not change (very much).  The current limit of max.iter =50 seems to work for most problems, but the Holzinger-Harmon 24 variable problem needs about 203 iterations to converge for a 5 factor solution.   

Principal axes may be used in cases when maximum likelihood solutions fail to converge. 

A problem in factor analysis is to find the best estimate of the original communalities.  Using the Squared Multiple Correlation (SMC) for each variable will underestimate the communalities, using 1s will over estimate.  By default, the SMC estimate is used.  If, however, a solution fails to be achieved, it is useful to try again using ones (SMC =FALSE).

The algorithm does not attempt to find the best (as defined by a maximum likelihood criterion) solution, but rather one that converges rapidly using successive eigen value decompositions.  The maximum likelihood criterion of fit and the associated chi square value are reported, and will be worse than that found using maximum likelihood procedures.

Although for items, it is typical to find factor scores by scoring the salient items (using, e.g.,\code{\link{score.items}} factor scores can be estimated by regression.  
}
\value{
  \item{values }{Eigen values of the final solution}
   \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}
  \item{rotation}{which rotation was requested?}
  \item{n.obs}{number of observations specified or found}
  \item{loadings}{An item by factor loading matrix of class ``loadings"  Suitable for use in other programs (e.g., GPA rotation or factor2cluster.}
  \item{fit}{How well does the factor model reproduce the correlation matrix. (See \code{\link{VSS}}, \code{\link{ICLUST}}, and \code{\link{principal}} for this fit statistic.}
  \item{fit.off}{how well are the off diagonal elements reproduced?}
  \item{dof}{Degrees of Freedom for this model. This is the number of observed correlations minus the number of independent parameters.  Let n=Number of items, nf = number of factors then
  \cr
  \eqn{dof = n * (n-1)/2 - n * nf + nf*(nf-1)/2}{dof = n * (n-1)/2 - n * nf + nf*(nf-1)/2}}
  \item{objective}{value of the function that is minimized by maximum likelihood procedures.  This is reported for comparison purposes and as a way to estimate chi square goodness of fit.  The objective function is 
  \cr
  \eqn{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}. }
  
  \item{STATISTIC}{If the number of observations is specified or found, this is a chi square based upon the objective function, f.  Using the formula from \code{\link{factanal}}(which seems to be Bartlett's test) :
 \cr
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } }
 
 \item{PVAL}{If n.obs > 0, then what is the probability of observing a chisquare this large or larger?}
 
 \item{Phi}{If oblique rotations (using oblimin from the GPArotation package or promax) are requested, what is the interfactor correlation.}
 \item{communality.iterations}{The history of the communality estimates.  Probably only useful for teaching what happens in the process of iterative fitting.}
 \item{residual}{If residuals are requested, this is the matrix of residual correlations after the factor model is applied.} 
 }
\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book.html} 
}
\author{ William Revelle }

\seealso{ \code{\link{principal}}, \code{\link{VSS}}, \code{\link{ICLUST}}}
\examples{
#using the Harman 24 mental tests, compare a principal factor with a principal components solution
pc <- principal(Harman74.cor$cov,4,rotate="varimax")
pa <- factor.pa(Harman74.cor$cov,4,rotate="varimax")
round(factor.congruence(pc,pa),2)

#then compare with a maximum likelihood solution using factanal
mle <- factanal(x,4,covmat=Harman74.cor$cov)
round(factor.congruence(mle,pa),2)
#note that the order of factors and the sign of some of factors differ 

#finally, compare the unrotated factor and pca solutions
pc1 <- principal(Harman74.cor$cov,4,rotate="none")
pa1 <- factor.pa(Harman74.cor$cov,4,rotate="none")
round(factor.congruence(pc1,pa1),2)
#note that the order of factors and the sign of some of factors differ 
}

\keyword{ multivariate }
\keyword{ models}% __ONLY ONE__ keyword per line
