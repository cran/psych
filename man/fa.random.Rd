\name{fa.random}
\alias{fa.random}

\title{A first approximation to Random Effects Exploratory Factor Analysis}

\description{Inspired, in part, by the wprifm function in the profileR package, fa.random removes between subject differences in mean level and then does a normal exploratory factor analysis of the ipsatized data.  Functionally, this removes a general factor of the data before factoring. To prevent non-positive definiteness of the residual data matrix, a very small amount of random noise is added to each variable. This is just a call to fa after removing the between subjects effect. Read the help file for \code{\link{fa}} for a detailed explanation of all of the input parameters and the output objects. 
}
\usage{
fa.random(data, nfactors = 1, fix = TRUE, n.obs = NA, n.iter = 1, rotate = "oblimin",
 scores = "regression", residuals = FALSE, SMC = TRUE, covar = FALSE, missing = FALSE,  
 impute = "median", min.err = 0.001, max.iter = 50, symmetric = TRUE, warnings = TRUE,
  fm = "minres", alpha = 0.1, p = 0.05, oblique.scores = FALSE, np.obs = NULL, 
  use = "pairwise", cor = "cor", weight = NULL, ...)
}
\arguments{
  \item{data}{A raw data matrix (or data.frame)}
  \item{nfactors}{ Number of factors to extract, default is 1 }
  \item{fix}{If TRUE, then a small amount of random error is added to each observed variable to keep the matrix positive semi-definite.  If FALSE, then this is not done but because the matrix is non-positive semi-definite it will need to be smoothed when finding the scores and the various statistics.}
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals. Ignored.}
  \item{np.obs}{The pairwise number of observations. Used if using a correlation matrix and asking for a minchi solution.}
   \item{rotate}{"none", "varimax", "quartimax",  "bentlerT", "equamax", "varimin", "geominT" and "bifactor" are orthogonal rotations.  "Promax", "promax", "oblimin", "simplimax", "bentlerQ,  "geominQ" and "biquartimin" and "cluster" are possible oblique transformations of the solution.  The default is to do a oblimin transformation, although  versions prior to 2009 defaulted to varimax. SPSS seems to do a Kaiser normalization before doing Promax, this is done here by the call to "promax" which does the normalization before calling Promax in GPArotation.}
 \item{n.iter}{Number of bootstrap interations to do in fa or fa.poly}
  \item{residuals}{Should the residual matrix be shown }
  \item{scores}{the default="regression" finds factor scores using regression.  Alternatives for  estimating factor scores include simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms ( \code{\link{factor.scores}}). Although scores="tenBerge" is probably preferred for most solutions, it will lead to problems with some improper correlation matrices.  }
  \item{SMC}{Use squared multiple correlations (SMC=TRUE) or use 1  as initial communality estimate. Try using 1 if imaginary eigen values are reported. If SMC is a vector of length the number of variables, then these values are used as starting values in the case of fm='pa'. }
  \item{covar}{if covar is TRUE, factor the covariance matrix, otherwise factor the correlation matrix}
  \item{missing}{if scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean}
  \item{impute}{"median" or "mean" values are used to replace missing values}
  \item{min.err}{Iterate until the change in communalities is less than min.err}
  \item{max.iter}{Maximum number of iterations for convergence }
  \item{symmetric}{symmetric=TRUE forces symmetry by just looking at the lower off diagonal values}
  \item{warnings}{warnings=TRUE => warn if number of factors is too many }
  \item{fm}{Factoring method  fm="minres" will do a minimum residual as will fm="uls".  Both of these use a first derivative.  fm="ols" differs very slightly from "minres" in that it minimizes the entire residual matrix using an OLS procedure but uses the empirical first derivative.  This will be slower.  fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample size weighted chi square when treating pairwise correlations with different number of subjects per pair. fm ="minrank" will do a minimum rank factor analysis. "old.min" will do minimal residual the way it was done prior to April, 2017 (see discussion below).}
  \item{alpha}{alpha level for the confidence intervals for RMSEA}
  \item{p}{if doing iterations to find confidence intervals, what probability values should be found for the confidence intervals}
  \item{oblique.scores}{When factor scores are found, should they be based on the structure matrix (default) or the pattern matrix (oblique.scores=TRUE).  } 
  \item{weight}{If not NULL, a vector of length n.obs that contains weights for each observation. The NULL case is equivalent to all cases being weighted 1.}
  \item{use}{How to treat missing data, use="pairwise" is the default".  See cor for other options.}
  \item{cor}{How to find the correlations: "cor" is Pearson", "cov" is covariance, 
  "tet" is tetrachoric, "poly" is polychoric, "mixed" uses mixed cor for a mixture of tetrachorics, polychorics, Pearsons, biserials, and polyserials, Yuleb is Yulebonett, Yuleq and YuleY are the obvious Yule coefficients as appropriate}

  \item{...}{additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax}
  
 } 

\details{This function is inspired by the wprifm function in the profileR package and the citation there to a paper by Davison, Kim and Close (2009).  The basic logic is to extract a means vector from each subject and then to analyze the resulting ipsatized data matrix.  This can be seen as removing acquiecence in the case of personality items, or the general factor, in the case of ability items.  Factors composed of items that are all keyed the same way (e.g., Neuroticism in the \code{\link{bfi}} data set) will be most affected by this technique. 

The output is identical to the normal \code{link{fa}} output with the addition of two objects:  subject and within.r.  The subject object is just the vector of the mean score for each subject on all the items. within.r is just the correlation of each item with those scores.  
}
\value{
\item{subject}{A vector of the  mean score on all items for each subject}
 \item{within.r}{The correlation of each item with the subject vector}
  \item{values }{Eigen values of the common factor solution}
  \item{e.values}{Eigen values of the original matrix}
  \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}
  \item{communalities}{If using minrank factor analysis, these are the communalities reflecting the total amount of common variance.  They will exceed the communality (above) which is the model estimated common variance. }
  \item{rotation}{which rotation was requested?}
  \item{n.obs}{number of observations specified or found}
  \item{loadings}{An item by factor (pattern) loading matrix of class ``loadings"  Suitable for use in other programs (e.g., GPA rotation or factor2cluster. To show these by sorted order, use \code{\link{print.psych}} with sort=TRUE}
  \item{complexity}{Hoffman's index of complexity for each item.  This is just \eqn{\frac{(\Sigma a_i^2)^2}{\Sigma a_i^4}}{{(\Sigma a_i^2)^2}/{\Sigma a_i^4}} where a_i is the factor loading on the ith factor. From Hofmann (1978), MBR. See also  Pettersson and Turkheimer (2010).}
  \item{Structure}{An item by factor structure matrix of class ``loadings". This is just the loadings (pattern) matrix times the factor intercorrelation matrix.}
  \item{fit}{How well does the factor model reproduce the correlation matrix. This is just \eqn{\frac{\Sigma r_{ij}^2 - \Sigma r^{*2}_{ij} }{\Sigma r_{ij}^2}
}{(sum(r^2ij - sum(r*^2ij))/sum(r^2ij}  (See \code{\link{VSS}}, \code{\link{ICLUST}}, and \code{\link{principal}} for this fit statistic.}
  \item{fit.off}{how well are the off diagonal elements reproduced?}
  \item{dof}{Degrees of Freedom for this model. This is the number of observed correlations minus the number of independent parameters.  Let n=Number of items, nf = number of factors then
  \cr
  \eqn{dof = n * (n-1)/2 - n * nf + nf*(nf-1)/2}{dof = n * (n-1)/2 - n * nf + nf*(nf-1)/2}}
  \item{objective}{Value of the function that is minimized by a maximum likelihood procedures.  This is reported for comparison purposes and as a way to estimate chi square goodness of fit.  The objective function is 
  \cr
  \eqn{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}. When using MLE, this function is minimized.  When using OLS (minres), although we are not minimizing this function directly, we can still calculate it in order to compare the solution to a MLE fit. }
  
  \item{STATISTIC}{If the number of observations is specified or found, this is a chi square based upon the objective function, f (see above). Using the formula from \code{\link{factanal}}(which seems to be Bartlett's test) :
 \cr
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } }
 
 \item{PVAL}{If n.obs > 0, then what is the probability of observing a chisquare this large or larger?}
 
 \item{Phi}{If oblique rotations (e.g,m using oblimin from the GPArotation package or promax) are requested, what is the interfactor correlation?}
 \item{communality.iterations}{The history of the communality estimates (For principal axis only.) Probably only useful for teaching what happens in the process of iterative fitting.}
 
 \item{residual}{The matrix of residual correlations after the factor model is applied. To display it conveniently, use the \code{\link{residuals}} command. }

\item{chi}{When normal theory fails (e.g., in the case of non-positive definite matrices), it useful to examine the empirically derived \eqn{\chi^2}{chi^2} based upon the sum of the squared residuals * N.  This will differ slightly from the MLE estimate which is based upon the fitting function rather than the actual residuals.}
  
 \item{rms}{This is the sum of the squared (off diagonal residuals) divided by the degrees of freedom.  Comparable to an RMSEA which, because it is based upon  \eqn{\chi^2}{chi^2}, requires the number of observations to be specified.  The rms is an empirical value while the RMSEA is based upon normal theory and the non-central \eqn{\chi^2}{chi^2} distribution. That is to say, if the residuals are particularly non-normal, the rms value and the associated  \eqn{\chi^2}{chi^2} and RMSEA can differ substantially. } 
 
 \item{crms}{rms adjusted for degrees of freedom}
 
 \item{RMSEA}{The Root Mean Square Error of Approximation is based upon the non-central 
 \eqn{\chi^2}{chi^2} distribution and the \eqn{\chi^2}{chi^2} estimate found from the MLE fitting function.  With normal theory data, this is fine.  But when the residuals are not distributed according to a noncentral \eqn{\chi^2}{chi^2}, this can give very strange values.  (And thus the confidence intervals can not be calculated.) The RMSEA is a conventional index of goodness (badness) of fit but it is also useful to examine the actual rms values.  }   
 
 
\item{TLI}{The Tucker Lewis Index of factoring reliability which is also known as the non-normed fit index.  }

\item{BIC}{Based upon \eqn{\chi^2}{chi^2} with the assumption of normal theory and using the \eqn{\chi^2}{chi^2} found using the objective function defined above. This is just \eqn{\chi^2 - 2 df}{chi^2 - 2 df}}


\item{eBIC}{When normal theory fails (e.g., in the case of non-positive definite matrices), it useful to examine the empirically derived eBIC based upon the empirical \eqn{\chi^2}{chi^2} - 2 df. }

\item{R2}{The multiple R square between the factors and factor score estimates, if they were to be found. (From Grice, 2001).  Derived from R2 is is the minimum correlation between any two factor estimates = 2R2-1. }

\item{r.scores}{The correlations of the factor score estimates using the specified model, if they were to be found.  Comparing these correlations with that of the scores themselves will show, if an alternative estimate of factor scores is used (e.g., the tenBerge method), the problem of factor indeterminacy.  For these correlations will not necessarily be the same.  }

\item{weights}{The beta weights to find the factor score estimates.  These are also used by the \code{\link{predict.psych}} function to find predicted factor scores for new cases.  These weights will depend upon the scoring method requested.  }
\item{scores}{The factor scores as requested.  Note that these scores reflect the choice of the way scores should be estimated (see scores in the input).  That is, simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms (see \code{\link{factor.scores}}).  The correlation between factor score estimates (r.scores) is based upon using the regression/Thurstone approach.  The actual correlation between scores will reflect the rotation algorithm chosen and may be found by correlating those scores. Although the scores are found by multiplying the standarized data by the weights matrix, this will not result in standard scores if using regression. }

\item{valid}{The validity coffiecient of course coded (unit weighted) factor score estimates (From Grice, 2001)}

\item{score.cor}{The correlation matrix of course coded (unit weighted) factor score estimates, if they were to be found, based upon the loadings matrix rather than the weights matrix.  }

\item{rot.mat}{The rotation matrix as returned from GPArotation.}
 }

\references{
 Davison, Mark L. and  Kim, Se-Kang and  Close, Catherine (2009)  Factor Analytic Modeling of Within Person Variation in Score Profiles.  Multivariate Behavioral Research (44(5) 668-687.

}
\author{
William Revelle
}
\note{An interesting, but not necessarily good, idea. 
To see what this does if there is a general factor, consider the unrotated solutions to the ability data set. In particular, compare the first factor loading with its congruence to the ipsatized solution means vector correlated with the items (the within.r object).  }

\seealso{
 \code{\link{fa}}
}
\examples{
fa.ab <- fa(ability,4,rotate="none")  #normal factor analysis
fa.ab.ip <- fa.random(ability,3,rotate="none") 
fa.congruence(list(fa.ab,fa.ab.ip,fa.ab.ip$within.r))


  }

\keyword{ multivariate }
\keyword{ models}