\name{principal}
\alias{principal}
\title{ Principal components analysis }
\description{Does an eigen value decomposition and returns eigen values, loadings, and degree of fit for a specified number of components.  Basically it is just  doing a principal components for n principal components.  Can show the residual correlations as well. The quality of reduction in the squared correlations is reported by comparing residual correlations to original correlations. Unlike princomp, this returns a subset of just the best nfactors. The eigen vectors are rescaled by the sqrt of the eigen values to produce the component loadings more typical in factor analysis.
}
\usage{
principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, scores=FALSE,missing=FALSE,impute="median")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a correlation matrix.  If a raw data matrix is used, the correlations will be found using pairwise deletions for missing values.}
  \item{nfactors}{Number of components to extract }
  \item{residuals}{ FALSE, do not show residuals, TRUE, report residuals }
  \item{rotate}{"none", "varimax", "quatimax", "promax", "oblimin", "simplimax", and "cluster" are possible rotations/transformations of the solution.}
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics.}
  \item{scores}{If TRUE, find component scores}
  \item{missing}{if scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean}
  \item{impute}{"median" or "mean" values are used to replace missing values} 
}
\details{Useful for those cases where the correlation matrix is improper (perhaps because of SAPA techniques).

There are a number of data reduction techniques including principal components and factor analysis.  Both PC and FA attempt to approximate a given correlation or covariance matrix of rank n with matrix of lower rank (p).  \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2} where k is much less than n.  For principal components, the item uniqueness is assumed to be zero and all elements of the correlation matrix are fitted. That is, \eqn{_nR_n \approx _{n}F_{kk}F_n'}{nRn = nFk kFn' }  The primary empirical difference between a components versus a factor model is the treatment of the variances for each item.  Philosophically, components are weighted composites of observed variables while in the factor model, variables are weighted composites of the factors.

For a n x n correlation matrix, the n principal components completely reproduce the correlation matrix.  However, if just the first k principal components are extracted, this is the best k dimensional approximation of the matrix.

It is important to recognize that rotated principal components are not principal components (the axes associated with the eigen value decomposition) but are merely components.  To point this out, unrotated principal components are labelled as PCi, while rotated PCs are now labeled as RCi (for rotated components) and obliquely transformed components as TCi (for transformed components). (Thanks to Ulrike Gromping for this suggestion.)

Rotations and transformations are either part of psych (Promax and cluster), of base R (varimax), or of GPArotation (simplimax, quartimax, oblimin).

Some of the statistics reported are more appropriate for (maximum likelihood) factor analysis rather than principal components analysis, and are reported to allow comparisons with these other models. 

Although for items, it is typical to find component scores by scoring the salient items (using, e.g., \code{\link{score.items}}) component scores can be estimated either by regression (which is done) or just multiplying the data by  the components.  The regression approach is done  to be parallel with the factor analysis function \code{\link{fa}}.  The regression weights are found from the inverse of the correlation matrix times the component loadings.  This has the result that the component scores are standard scores (mean=0, sd = 1) of the standardized input.  A comparison to the scores from \code{\link{princomp}} shows this difference.  princomp does not, by default, standardize the data matrix, nor are the components themselves standardized.    
}
\value{
  
  \item{values}{Eigen Values of all components -- useful for a scree plot}
  \item{rotation}{which rotation was requested?}
  \item{n.obs}{number of observations specified or found}
   \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}

  \item{loadings }{A standard loading matrix of class ``loadings"}
  \item{fit }{Fit of the model to the correlation matrix }
  \item{fit.off}{how well are the off diagonal elements reproduced?}
  \item{residual }{Residual matrix -- if requested}
  \item{dof}{Degrees of Freedom for this model. This is the number of observed correlations minus the number of independent parameters (number of items * number of factors - nf*(nf-1)/2.   That is, dof = niI * (ni-1)/2 - ni * nf + nf*(nf-1)/2.}
  \item{objective}{value of the function that is minimized by maximum likelihood procedures.  This is reported for comparison purposes and as a way to estimate chi square goodness of fit.  The objective function is 
  \cr
  \eqn{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}.   Because components do not minimize the off diagonal, this fit will be not as good as for factor analysis. }
  
  \item{STATISTIC}{If the number of observations is specified or found, this is a chi square based upon the objective function, f.  Using the formula from \code{\link{factanal}}:
 \cr
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } }
 
 \item{PVAL}{If n.obs > 0, then what is the probability of observing a chisquare this large or larger?}
 \item{phi}{If oblique rotations (using oblimin from the GPArotation package) are requested, what is the interfactor correlation.}
 \item{scores}{If scores=TRUE, then estimates of the factor scores are reported }
 \item{weights}{The beta weights to find the principal components from the data}
 \item{R2}{The multiple R square between the factors and factor score estimates, if they were to be found. (From Grice, 2001)  For components, these are of course 1.0.}
\item{valid}{The correlations of the component score estimates with the components, if they were to be found and unit weights were used. (So called course coding).}
}

\author{ William Revelle}
\references{

Grice, James W.  (2001), Computing and evaluating factor scores. Psychological Methods, 6, 430-450

Revelle, W. An introduction to psychometric theory with applications in R (in prep) Springer. Draft chapters available at  \url{http://personality-project.org/r/book/}
}

\seealso{\code{\link{VSS}} (to test for the number of components or factors to extract), \code{\link{VSS.scree}} and \code{\link{fa.parallel}} to show a scree plot and compare it with random resamplings of the data), \code{\link{factor2cluster}} (for course coding keys), \code{\link{fa}} (for factor analysis), \code{\link{factor.congruence}} (to compare solutions)}
\examples{
#Four principal components of the Harmon 24 variable problem
#compare to a four factor principal axes solution using factor.congruence
pc <- principal(Harman74.cor$cov,4,rotate="varimax")
mr <- fa(Harman74.cor$cov,4,rotate="varimax")  #minres factor analysis
pa <- fa(Harman74.cor$cov,4,rotate="varimax",fm="pa")  # principal axis factor analysis
round(factor.congruence(list(pc,mr,pa)),2)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
