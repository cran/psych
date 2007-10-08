\name{principal}
\alias{principal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Principal components analysis }
\description{Does an eigen value decomposition and returns eigen values, loadings, and degree of fit for a specified number of components.  Basically it is just  doing a principal components for n principal components.  Can show the residual correlations as well. The quality of reduction in the squared correlations is reported by comparing residual correlations to original correlations. Unlike princomp, this returns a subset of just the best nfactors. The eigen vectors are rescaled by the sqrt of the eigen values to produce the component loadings more typical in factor analysis.
}
\usage{
principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NULL, digits=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a correlation matrix.  If a raw data matrix is used, the correlations will be found using pairwise deletions for missing values.}
  \item{nfactors}{Number of components to extract }
  \item{residuals}{ FALSE, do not show residuals, TRUE, report residuals }
  \item{rotate} {"none", "varimax", "promax" or "oblimin" are possible rotations of the solution.}
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics.}
  \item{digits }{digits =2.    Accuracy of answers as well as display}
}
\details{Useful for those cases where the correlation matrix is improper (perhaps because of SAPA techniques).

There are a number of data reduction techniques including principal components and factor analysis.  Both PC and FA attempt to approximate a given correlation or covariance matrix of rank n with matrix of lower rank (p).  \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2} where k is much less than n.  For principal components, the item uniqueness is assumed to be zero and all elements of the correlation matrix are fitted. That is, \eqn{_nR_n \approx _{n}F_{kk}F_n'}{nRn = nFk kFn' }  The primary empirical difference between a components versus a factor model is the treatment of the variances for each item.  Philosophically, components are weighted composites of observed variables while in the factor model, variables are weighted composites of the factors.

For a n x n correlation matrix, the n principal components completely reproduce the correlation matrix.  However, if just the first k principal components are extracted, this is the best k dimensional approximation of the matrix.

Some of the statistics reported are more appropriate for maximum likelihood factor analysis rather than principal components analysis, and are reported to allow comparisons with these other models. 
}
\value{
  
  \item{values }{Eigen Values of all components -- useful for a scree plot}
  \item{rotation}{which rotation was requested?}
  \item{n.obs}{number of observations specified or found}
   \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}

  \item{loadings }{A standard loading matrix of class ``loadings"}
  \item{fit }{Fit of the model to the correlation matrix }
   \item{fit.off}{how well are the off diagonal elements reproduced?}
  \item{residual }{Residual matrix -- if requested}
   \item{communality}{The history of the communality estimates.  Probably only useful for teaching what happens in the process of iterative fitting.}
  \item{dof}{Degrees of Freedom for this model. This is the number of observed correlations minus the number of independent parameters (number of items * number of factors - nf*(nf-1)/2.   That is, dof = niI * (ni-1)/2 - ni * nf + nf*(nf-1)/2.}
  \item{objective}{value of the function that is minimized by maximum likelihood procedures.  This is reported for comparison purposes and as a way to estimate chi square goodness of fit.  The objective function is 
  \cr
  \eqn{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}. }
  
  \item{STATISTIC}{If the number of observations is specified or found, this is a chi square based upon the objective function, f.  Using the formula from \code{\link{factanal}}:
 \cr
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } }
 
 \item{PVAL}{If n.obs > 0, then what is the probability of observing a chisquare this large or larger?}
 \item{phi}{If oblique rotations (using oblimin from the GPArotation package) are requested, what is the interfactor correlation.}
}

\author{ William Revelle}

\seealso{\code{\link{VSS}},\code{\link{factor2cluster}},\code{\link{factor.pa}}, \code{\link{factor.congruence}}}
\examples{
#Four principal components of the Harmon 24 variable problem
#compare to a four factor principal axes solution using factor.congruence
pc <- principal(Harman74.cor$cov,4,rotate=TRUE)
pa <- factor.pa(Harman74.cor$cov,4,rotate=TRUE)
round(factor.congruence(pc,pa),2)


}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
