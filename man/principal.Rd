\name{principal}
\alias{principal}
\alias{pca}
\title{ Principal components analysis (PCA)}
\description{Does an eigen value decomposition and returns eigen values, loadings, and degree of fit for a specified number of components.  Basically it is just  doing a principal components analysis (PCA) for n principal components of either a correlation or covariance matrix.  Can show the residual correlations as well. The quality of reduction in the squared correlations is reported by comparing residual correlations to original correlations. Unlike princomp, this returns a subset of just the best nfactors. The eigen vectors are rescaled by the sqrt of the eigen values to produce the component loadings more typical in factor analysis.
}
\usage{
principal(r, nfactors = 1, residuals = FALSE,rotate="varimax",n.obs=NA, covar=FALSE,
 scores=TRUE,missing=FALSE,impute="median",oblique.scores=TRUE,method="regression",
 use ="pairwise",cor="cor",correct=.5,weight=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a correlation matrix.  If a raw data matrix is used, the correlations will be found using pairwise deletions for missing values.}
  \item{nfactors}{Number of components to extract }
  \item{residuals}{ FALSE, do not show residuals, TRUE, report residuals }
  \item{rotate}{"none", "varimax", "quartimax", "promax", "oblimin", "simplimax", and "cluster" are possible rotations/transformations of the solution. See \code{\link{fa}} for all rotations avaiable.}
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics.}
  \item{covar}{If false, find the correlation matrix from the raw data or convert to a correlation matrix if given a square matrix as input.}
  \item{scores}{If TRUE, find component scores}
  \item{missing}{if scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean}
  \item{impute}{"median" or "mean" values are used to replace missing values} 
  \item{oblique.scores}{If TRUE (default), then the component scores are based upon the structure matrix.  If FALSE, upon the pattern matrix.}
  \item{method}{Which way of finding component scores should be used. The default is "regression"}
   \item{weight}{If not NULL, a vector of length n.obs that contains weights for each observation. The NULL case is equivalent to all cases being weighted 1.}
  \item{use}{How to treat missing data, use="pairwise" is the default".  See cor for other options.}
\item{cor}{How to find the correlations: "cor" is Pearson", "cov" is covariance, 
  "tet" is tetrachoric, "poly" is polychoric, "mixed" uses mixedCor for a mixture of tetrachorics, polychorics, Pearsons, biserials, and polyserials, Yuleb is Yulebonett, Yuleq and YuleY are the obvious Yule coefficients as appropriate}
\item{correct}{When doing tetrachoric, polycoric, or mixed cor, how should we treat empty cells.  (See the discussion in the help for tetrachoric.)} 
  \item{...}{other parameters to pass to functions such as factor.scores or the various rotation functions.  }
  }
\details{Useful for those cases where the correlation matrix is improper (perhaps because of SAPA techniques).

There are a number of data reduction techniques including principal components analysis (PCA) and factor analysis (EFA).  Both PC and FA attempt to approximate a given correlation or covariance matrix of rank n with matrix of lower rank (p).  \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2} where k is much less than n.  For principal components, the item uniqueness is assumed to be zero and all elements of the correlation or covariance matrix are fitted. That is, \eqn{_nR_n \approx _{n}F_{kk}F_n'}{nRn = nFk kFn' }  The primary empirical difference between a components versus a factor model is the treatment of the variances for each item.  Philosophically, components are weighted composites of observed variables while in the factor model, variables are weighted composites of the factors. As the number of items increases, the difference between the two models gets smaller.  Factor loadings are the asymptotic component loadings as the number of items gets larger. 

For a n x n correlation matrix, the n principal components completely reproduce the correlation matrix.  However, if just the first k principal components are extracted, this is the best k dimensional approximation of the matrix.

It is important to recognize that rotated principal components are not principal components (the axes associated with the eigen value decomposition) but are merely components.  To point this out, unrotated principal components are labelled as PCi, while rotated PCs are now labeled as RCi (for rotated components) and obliquely transformed components as TCi (for transformed components). (Thanks to Ulrike Gromping for this suggestion.)

Rotations and transformations are either part of psych (Promax and cluster), of base R (varimax), or of GPArotation (simplimax, quartimax, oblimin, etc.).  

Of the various rotation/transformation options, varimax, Varimax, quartimax, bentlerT, geominT, and bifactor do orthogonal rotations. Promax  transforms obliquely with a target matix equal to the varimax solution. oblimin, quartimin, simplimax,  bentlerQ,  geominQ and biquartimin are oblique transformations. Most of these are just calls to  the GPArotation package. The ``cluster'' option does a targeted rotation to a structure defined by the cluster representation of a varimax solution.  With the optional "keys" parameter, the "target" option will rotate to a target supplied as a keys matrix. (See \code{\link{target.rot}}.)

The rotation matrix  (rot.mat) is returned from all of these options. This is the inverse of the Th (theta?) object returned by the GPArotation package.  The correlations of the factors may be found by \eqn{\Phi = \theta' \theta}{Phi = Th' Th}

Some of the statistics reported are more appropriate for (maximum likelihood) factor analysis rather than principal components analysis, and are reported to allow comparisons with these other models. 

Although for items, it is typical to find component scores by scoring the salient items (using, e.g., \code{\link{scoreItems}}) component scores are found  by regression where the regression weights are \eqn{R^{-1} \lambda}{R^(-1) lambda} where \eqn{\lambda}{lambda} is the matrix of component loadings.   The regression approach is done  to be parallel with the factor analysis function \code{\link{fa}}.  The regression weights are found from the inverse of the correlation matrix times the component loadings.   This has the result that the component scores are standard scores (mean=0, sd = 1) of the standardized input.  A comparison to the scores from \code{\link{princomp}} shows this difference.  princomp does not, by default, standardize the data matrix, nor are the components themselves standardized.  The regression weights are found from the Structure matrix, not the Pattern matrix. If the scores are found with the covar option = TRUE, then the scores are not standardized but are just mean centered.  

Jolliffe (2002) discusses why the interpretation of rotated components is complicated.   Rencher (1992) discourages the use of rotated components. The approach used here is consistent with the factor analytic tradition.  The correlations of the items with the component scores closely matches (as it should) the component loadings (as reported in the structure matrix). 

The output from the print.psych function displays the component loadings (from the pattern matrix), the h2 (communalities) the u2 (the uniquenesses), com (the complexity of the component loadings for that variable (see below).  In the case of an orthogonal solution, h2 is merely the row sum of the squared component loadings. But for an oblique solution, it is the row sum of the (squared) orthogonal component loadings (remember, that rotations or transformations do not change the communality).  This information is returned (invisibly) from the print function as the object Vaccounted.
}
\value{
  
  \item{values}{Eigen Values of all components -- useful for a scree plot}
  \item{rotation}{which rotation was requested?}
  \item{n.obs}{number of observations specified or found}
   \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}
    \item{complexity}{Hoffman's index of complexity for each item.  This is just \eqn{\frac{(\Sigma a_i^2)^2}{\Sigma a_i^4}}{{(\Sigma a_i^2)^2}/{\Sigma a_i^4}} where a_i is the factor loading on the ith factor. From Hofmann (1978), MBR. See also  Pettersson and Turkheimer (2010).}

  \item{loadings }{A standard loading matrix of class ``loadings"}
  \item{fit }{Fit of the model to the correlation matrix }
  \item{fit.off}{how well are the off diagonal elements reproduced?}
  \item{residual }{Residual matrix -- if requested}
  \item{dof}{Degrees of Freedom for this model. This is the number of observed correlations minus the number of independent parameters (number of items * number of factors - nf*(nf-1)/2.   That is, dof = niI * (ni-1)/2 - ni * nf + nf*(nf-1)/2.}
  \item{objective}{value of the function that is minimized by maximum likelihood procedures.  This is reported for comparison purposes and as a way to estimate chi square goodness of fit.  The objective function is 
  \cr
  \eqn{f = (trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}.   Because components do not minimize the off diagonal, this fit will be not as good as for factor analysis. It is included merely for comparison purposes.}
  
  \item{STATISTIC}{If the number of observations is specified or found, this is a chi square based upon the objective function, f.  Using the formula from \code{\link{factanal}}:
 \cr
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } }
 
 \item{PVAL}{If n.obs > 0, then what is the probability of observing a chisquare this large or larger?}
 \item{Phi}{If oblique rotations (using oblimin from the GPArotation package) are requested, what is the interfactor correlation.}
 \item{scores}{If scores=TRUE, then estimates of the factor scores are reported }
 \item{weights}{The beta weights to find the principal components from the data}
 \item{R2}{The multiple R square between the factors and factor score estimates, if they were to be found. (From Grice, 2001)  For components, these are of course 1.0.}
\item{valid}{The correlations of the component score estimates with the components, if they were to be found and unit weights were used. (So called course coding).}
\item{r.scores}{The correlation of the component scores.  (Since components are just weighted linear sums of the items, this is the same as Phi).}
\item{rot.mat}{The rotation matrix used to produce the rotated component loadings. }
}


\note{By default, the accuracy of the varimax rotation function seems to be less than the Varimax function. This can be enhanced by specifying eps=1e-14 in the call to principal if using varimax rotation. Furthermore, note that Varimax by default does not apply the Kaiser normalization, but varimax does. Gottfried Helms compared these two rotations with those produced by SPSS and found identical values if using the appropriate options. (See the last two examples.)

The ability to use different kinds of correlations was added in version 1.9.12.31 to be compatible with the options in fa.

}
\author{ William Revelle}
\references{

Grice, James W.  (2001), Computing and evaluating factor scores. Psychological Methods, 6, 430-450

Jolliffe, I. (2002)  Principal Component Analysis (2nd ed). Springer.

Rencher, A. C. (1992) Interpretation of Canonical Discriminant Functions, Canonical Variates, and Principal Components, the American Statistician, (46) 217-225.

Revelle, W. An introduction to psychometric theory with applications in R (in prep) Springer. Draft chapters available at  \url{https://personality-project.org/r/book/}
}

\seealso{\code{\link{VSS}} (to test for the number of components or factors to extract), \code{\link{VSS.scree}} and \code{\link{fa.parallel}} to show a scree plot and compare it with random resamplings of the data), \code{\link{factor2cluster}} (for course coding keys), \code{\link{fa}} (for factor analysis), \code{\link{factor.congruence}} (to compare solutions), \code{\link{predict.psych}} to find factor/component scores for a new data set based upon the weights from an original data set. }
\examples{
#Four principal components of the Harman 24 variable problem
#compare to a four factor principal axes solution using factor.congruence
pc <- principal(Harman74.cor$cov,4,rotate="varimax")
mr <- fa(Harman74.cor$cov,4,rotate="varimax")  #minres factor analysis
pa <- fa(Harman74.cor$cov,4,rotate="varimax",fm="pa")  # principal axis factor analysis
round(factor.congruence(list(pc,mr,pa)),2)

pc2 <- principal(Harman.5,2,rotate="varimax")
pc2
round(cor(Harman.5,pc2$scores),2)  #compare these correlations to the loadings 
#now do it for unstandardized scores, and transform obliquely
pc2o <- principal(Harman.5,2,rotate="promax",covar=TRUE)
pc2o
round(cov(Harman.5,pc2o$scores),2) 
pc2o$Structure    #this matches the covariances with the scores
biplot(pc2,main="Biplot of the Harman.5 socio-economic variables",labels=paste0(1:12))

#For comparison with SPSS  (contributed by Gottfried Helms)
pc2v <- principal(iris[1:4],2,rotate="varimax",normalize=FALSE,eps=1e-14)
print(pc2v,digits=7)
pc2V <- principal(iris[1:4],2,rotate="Varimax",eps=1e-7)
p <- print(pc2V,digits=7)
round(p$Vaccounted,2)   # the amount of variance accounted for is returned as an object of print
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
