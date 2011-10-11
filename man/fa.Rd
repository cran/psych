\name{fa}
\alias{fa}
\alias{factor.pa} 
\alias{factor.minres} 
\alias{factor.wls}  
\alias{fac}
\alias{fa.poly}



\title{MinRes (minimum residual) Factor analysis as well as Factor Analysis by Principal Axis, Weighted Least Squares or Maximum Likelihood }
\description{Among the many ways to do latent variable factor analysis, one of the better is to use Ordinary Least Squares to find the minimum residual (minres) solution. This produces solutions very similar to maximum likelihood even for badly behaved matrices. A variation on minres is to do weighted least squares. Perhaps the most conventional technique is principal axes.  An eigen value decomposition of a correlation matrix is done and then the communalities for each variable are estimated by the first n factors. These communalities are entered onto the diagonal and the procedure is repeated until the sum(diag(r)) does not vary.   Yet another estimate procedure is maximum likelihood. For well behaved matrices, maximum likelihood factor analysis (either in the fa or in the factanal fuction) is probably preferred.  Bootstrapped confidence intervals the loadings and interfactor correlations are found by fa  in n.iter > 1 and for dichotomous or polytomous items, by fa.poly.
}
\usage{
fa(r,nfactors=1,n.obs = NA,n.iter=1, rotate="oblimin", scores="regression", residuals=FALSE, SMC=TRUE, covar=FALSE,missing=FALSE,impute="median",min.err = 0.001,  max.iter = 50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1,p=.05,oblique.scores=TRUE,...)

fac(r,nfactors=1,n.obs = NA, rotate="oblimin", scores="tenBerge", residuals=FALSE, SMC=TRUE, covar=FALSE,missing=FALSE,impute="median",min.err = 0.001,  max.iter = 50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1,oblique.scores=TRUE,...)

fa.poly(x,nfactors=1,n.obs = NA, n.iter=1, rotate="oblimin", SMC=TRUE,  missing=FALSE, impute="median", min.err = .001, max.iter=50, symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05, oblique.scores=TRUE,...)

factor.minres(r, nfactors=1, residuals = FALSE, rotate = "varimax",n.obs = NA,
scores = FALSE,SMC=TRUE, missing=FALSE,impute="median",min.err = 0.001, digits = 2, max.iter = 50,symmetric=TRUE,warnings=TRUE,fm="minres")    #deprecated

factor.wls(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NA,
scores=FALSE,SMC=TRUE,missing=FALSE,impute="median", min.err = .001,digits=2,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="wls")  #deprecated
}

\arguments{
  \item{r}{A correlation matrix or a raw data matrix. If raw data, the correlation matrix will be found using pairwise deletion. }
  \item{x}{For fa.poly.ci, only raw data may be used}
  \item{nfactors}{ Number of factors to extract, default is 1 }
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.}
   \item{rotate}{"none", "varimax", "quartimax",  "bentlerT",  "geominT" and "bifactor" are orthogonal rotations.  "promax", "oblimin", "simplimax", "bentlerQ,  "geominQ" and "biquartimin" and "cluster" are possible rotations or transformations of the solution.  The default is to do a oblimin transformation, although  versions prior to 2009 defaulted to varimax.}
 \item{n.iter}{Number of bootstrap interations to do in fa or fa.poly}
  \item{residuals}{Should the residual matrix be shown }
  \item{scores}{the default="regression" finds factor scores using regression.  Alternatives for  estimating factor scores include simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms (see factor.scores). Although scores="tenBerge" is probably preferred for most solutions, it will lead to problems with some improper correlation matrices.  }
  \item{SMC}{Use squared multiple correlations (SMC=TRUE) or use 1  as initial communality estimate. Try using 1 if imaginary eigen values are reported. }
  \item{covar}{if covar is TRUE, factor the covariance matrix, otherwise factor the correlation matrix}
  \item{missing}{if scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean}
  \item{impute}{"median" or "mean" values are used to replace missing values}
  \item{min.err}{Iterate until the change in communalities is less than min.err}
  \item{digits}{ How many digits of output should be returned-- deprecated -- now specified in the print function}
  \item{max.iter}{Maximum number of iterations for convergence }
  \item{symmetric}{symmetric=TRUE forces symmetry by just looking at the lower off diagonal values}
  \item{warnings}{warnings=TRUE => warn if number of factors is too many }
  \item{fm}{factoring method  fm="minres" will do a minimum residual (OLS), fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis}
  \item{alpha}{alpha level for the confidence intervals for RMSEA}
  \item{p}{if doing iterations to find confidence intervals, what probability values should be found for the confidence intervals}
  \item{oblique.scores}{When factor scores are found, should they be based on the structure matrix (default) or the pattern matrix (which is what factanal seems to do). } 
  \item{...}{additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax}
  
 } 
 
\details{Factor analysis is an attempt to approximate a correlation or covariance matrix with one of lesser rank.  The basic model is that \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2} where k is much less than n. There are many ways to do factor analysis, and maximum likelihood procedures are probably the most preferred (see \code{\link{factanal}} ).  The existence of uniquenesses is what distinguishes factor analysis from principal components analysis (e.g., \code{\link{principal}}). If variables are thought to represent a ``true" or latent part then factor analysis provides an estimate of the correlations with the latent factor(s) representing the data.  If variables are thought to be measured without error, then principal components provides the most parsimonious description of the data.  

The fa function will do factor analyses using one of four different algorithms: minimum residual (minres), principal axes, weighted least squares, or maximum likelihood.

Principal axes factor analysis has a long history in exploratory analysis and is a straightforward procedure.  Successive eigen value decompositions are done on a correlation matrix with the diagonal replaced with  diag (FF') until sum(diag(FF')) does not change (very much).  The current limit of max.iter =50 seems to work for most problems, but the Holzinger-Harmon 24 variable problem needs about 203 iterations to converge for a 5 factor solution.  

Not all factor programs that do principal axes do iterative solutions.  The example from the SAS manual (Chapter 26) is such a case. To achieve that solution, it is necessary to specify that the max.iterations = 1.  Comparing that solution to an iterated one (the default) shows that iterations improve the solution.  In addition, fm="minres" or fm="mle" produces even better solutions for this example.  

Principal axes may be used in cases when maximum likelihood solutions fail to converge, although fm="minres" will also do that and tends to produce better (smaller residuals) solutions.

A problem in factor analysis is to find the best estimate of the original communalities.  Using the Squared Multiple Correlation (SMC) for each variable will underestimate the communalities, using 1s will over estimate.  By default, the SMC estimate is used.  In either case, iterative techniques will tend to converge on a stable solution. If, however, a solution fails to be achieved, it is useful to try again using ones (SMC =FALSE).

The iterated principal axes algorithm does not attempt to find the best (as defined by a maximum likelihood criterion) solution, but rather one that converges rapidly using successive eigen value decompositions.  The maximum likelihood criterion of fit and the associated chi square value are reported, and will be worse than that found using maximum likelihood procedures.

The minimum residual (minres) solution is an unweighted least squares solution that takes a slightly different approach.  It uses the \code{\link{optim}} function and adjusts the diagonal elements of the correlation matrix to mimimize the squared residual when the factor model is the eigen value decomposition of the reduced matrix.  MINRES and PA will both work when ML will not, for they can be used when the matrix is singular. At least on a number of test cases, the MINRES solution is slightly more similar to the ML solution than is the PA solution.  To a great extent, the minres and wls solutions follow ideas in the \code{\link{factanal}} function.  

The weighted least squares (wls) solution weights the residual matrix by 1/ diagonal of the inverse of the correlation matrix.  This has the effect of weighting items with low communalities more than those with high communalities. 

The generalized least squares (gls) solution weights the residual matrix by the inverse of the correlation matrix.  This has the effect of weighting those variables with low communalities even more than those with high communalities.

The maximum likelihood solution takes yet another approach and finds those communality values that minimize the chi square goodness of fit test.  The fm="ml" option provides a maximum likelihood solution following the procedures used in \code{\link{factanal}} but does not provide all the extra features of that function.

Test cases comparing the output to SPSS suggest that the PA algorithm matches what SPSS calls uls, and that the wls solutions are equivalent in their fits. The wls and gls solutions have slightly larger eigen values, but slightly worse fits of the off diagonal residuals than do the minres or maximum likelihood solutions.  Comparing the results to the examples in Harman (76), the PA solution with no iterations matches what Harman calls Principal Axes (as does SAS), while the iterated PA solution matches his minres solution.  The minres solution found in psych tends to have slightly smaller off diagonal residuals (as it should) than does the iterated PA solution.   

Although for items, it is typical to find factor scores by scoring the salient items (using, e.g., \code{\link{score.items}}) factor scores can be estimated by regression as well as several other means. There are multiple approaches that are possible (see Grice, 2001) and the one taken here was developed by tenBerge et al.  (see \code{\link{factor.scores}}. The alternative, which will match factanal is to find the scores using regression --  Thurstone's least squares regression where the weights are found by   \eqn{W = R^(-1)S}{W = inverse(R)S} where R is the correlation matrix of the variables ans S is the structure matrix.  Then, factor scores are just \eqn{Fs = X W}{Fs = X W}.

In the oblique case, the factor loadings are referred to as Pattern coefficients and are related to the Structure coefficients by \eqn{S = P \Phi}{S = P Phi} and thus \eqn{P = S \Phi^{-1}}{P = S Phi^-1}.  When estimating factor scores, \code{\link{fa}} and  \code{\link{factanal}} differ in that \code{\link{fa}} finds the factors from the Structure matrix and \code{\link{factanal}}  seems to find them from the Pattern matrix.  Thus, although in the orthogonal case, fa and factanal agree perfectly in their factor score estimates, they do not agree in the case of oblique factors.  Setting oblique.scores = FALSE will produce factor score estimate that match those of \code{\link{factanal}}.

It is sometimes useful to extend the factor solution to variables that were not factored.  This may be done using \code{\link{fa.extension}}. Factor extension is typically done in the case where some variables were not appropriate to factor, but factor loadings on the original factors are still desired.

For dichotomous items or polytomous items, it is recommended to analyze the \code{\link{tetrachoric}} or \code{\link{polychoric}} correlations rather than the Pearson correlations.  This is done automatically when using  \code{\link{irt.fa}} or \code{\link{fa.poly}}  functions.  In the first case,  the factor analysis results are reported in Item Response Theory (IRT) terms, although the original factor solution is returned in the results. In the later case, a typical factor loadings matrix is returned, but the tetrachoric/polychoric correlation matrix and item statitics are saved for reanalysis by \code{\link{irt.fa}}


Of the various rotation/transformation options, varimax, Varimax, quartimax, bentlerT, geominT, and bifactor do orthogonal rotations. Promax  transforms obliquely with a target matix equal to the varimax solution. oblimin, quartimin, simplimax,  bentlerQ,  geominQ and biquartimin are oblique transformations. Most of these are just calls to  the GPArotation package. The ``cluster'' option does a targeted rotation to a structure defined by the cluster representation of a varimax solution.  With the optional "keys" parameter, the "target" option will rotate to a target supplied as a keys matrix. (See \code{\link{target.rot}}.)

The "bifactor" rotation implements the Jennrich and Bentler (2011) bifactor rotation by calling the GPForth function in the GPArotation package and using two functions adapted from  the MatLab code of Jennrich and Bentler.

There are two varimax rotation functions.  One, Varimax, in the GPArotation package does not by default apply Kaiser normalization.  The other, varimax, in the stats package, does.  It appears that the two rotation functions produce slightly different results even when normalization is set. For consistency with the other rotation functions, Varimax is probably preferred.

When factor analyzing items with dichotomous or polytomous responses, the \code{\link{irt.fa}} function provides an Item Response Theory representation of the factor output. 

The function \code{\link{fa}} will repeat the analysis n.iter times on a bootstrapped sample of the data (if they exist) or of a simulated data set based upon the observed correlation matrix.  The mean estimate and standard deviation of the estimate are returned and will print the original factor analysis as well as the alpha level confidence intervals for the estimated coefficients.  The bootstrapped solutions are rotated towards the original solution using target.rot. The factor loadings are z-transformed, averaged and then back transformed. The default is to have n.iter =1 and thus not do bootstrapping.

\code{\link{fa.poly}} will find confidence intervals for a factor solution for dichotomous or polytomous items (set n.iter > 1 to do so).  Perhaps more useful is to find the Item Response Theory parameters equivalent to the factor loadings reported in fa.poly by using the \code{\link{irt.fa}} function.  

Some correlation matrices that arise from using pairwise deletion or from tetrachoric or polychoric matrices will not be proper.  That is, they will not be positive semi-definite (all eigen values >= 0).  The \code{\link{cor.smooth}} function will adjust correlation matrices (smooth them) by making all negative eigen values slightly greater than 0, rescaling the other eigen values to sum to the number of variables, and then recreating the correlation matrix.  See \code{\link{cor.smooth}} for an example of this problem using the \code{\link{burt}} data set.

}
\value{
  \item{values }{Eigen values of the common factor solution}
  \item{e.values}{Eigen values of the original matrix}
  \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}
  \item{rotation}{which rotation was requested?}
  \item{n.obs}{number of observations specified or found}
  \item{loadings}{An item by factor loading matrix of class ``loadings"  Suitable for use in other programs (e.g., GPA rotation or factor2cluster. To show these by sorted order, use \code{\link{print.psych}} with sort=TRUE}
  \item{fit}{How well does the factor model reproduce the correlation matrix. This is just \eqn{\frac{\Sigma r_{ij}^2 - \Sigma r^{*2}_{ij} }{\Sigma r_{ij}^2}
}{(sum(r^2ij - sum(r*^2ij))/sum(r^2ij}  (See \code{\link{VSS}}, \code{\link{ICLUST}}, and \code{\link{principal}} for this fit statistic.}
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
 \item{communality.iterations}{The history of the communality estimates (For principal axis only.) Probably only useful for teaching what happens in the process of iterative fitting.}
 \item{residual}{If residuals are requested, this is the matrix of residual correlations after the factor model is applied.} 
 
 \item{rms}{This is the sum of the squared (off diagonal residuals) divided by the degrees of freedom.  Comparable to an RMSEA which, because it is based upon  \eqn{\chi^2}{chi^2}, requires the number of observations to be specified.  The rms is an empirical value while the RMSEA is based upon normal theory and the non-central \eqn{\chi^2}{chi^2} distribution.}
 
\item{TLI}{The Tucker Lewis Index of factoring reliability which is also known as the non-normed fit index.  }

\item{BIC}{Based upon \eqn{\chi^2}{chi^2}}

\item{R2}{The multiple R square between the factors and factor score estimates, if they were to be found. (From Grice, 2001).  Derived from R2 is is the minimum correlation between any two factor estimates = 2R2-1. }

\item{r.scores}{The correlations of the factor score estimates, if they were to be found.}

\item{weights}{The beta weights to find the factor score estimates}

\item{valid}{The validity coffiecient of course coded (unit weighted) factor score estimates (From Grice, 2001)}

\item{score.cor}{The correlation matrix of course coded (unit weighted) factor score estimates, if they were to be found, based upon the loadings matrix rather than the weights matrix.  }

 }

\note{SPSS will sometimes use a Kaiser normalization before rotating.  This will lead to different solutions than reported here.  To get the Kaiser normalized loadings, use \code{\link{kaiser}}}

\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.

Grice, James W.  (2001), Computing and evaluating factor scores. Psychological Methods, 6, 430-450

Harman, Harry and Jones, Wayne (1966) Factor analysis by minimizing residuals (minres), Psychometrika, 31, 3, 351-368.


Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 
}
\author{ William Revelle }

\note{Thanks to Erich Studerus for some very helpful suggestions about various rotation and factor scoring algorithms, and to  Gumundur Arnkelsson  for suggestions about factor scores for singular matrices.

The fac function is the original fa function which is now called by fa repeatedly to get confidence intervals.  }

\seealso{ \code{\link{principal}},  \code{\link{irt.fa}}, \code{\link{VSS}}, \code{\link{ICLUST}}, \code{\link{predict.psych}}, \code{\link{fa.extension}}  }
\examples{
#using the Harman 24 mental tests, compare a principal factor with a principal components solution
pc <- principal(Harman74.cor$cov,4,rotate="varimax")
pa <- fa(Harman74.cor$cov,4,fm="pa" ,rotate="varimax")  #principal axis 
uls <- fa(Harman74.cor$cov,4,rotate="varimax")          #unweighted least squares is minres
wls <- fa(Harman74.cor$cov,4,fm="wls")       #weighted least squares

#to show the loadings sorted by absolute value
print(uls,sort=TRUE)

#then compare with a maximum likelihood solution using factanal
mle <- factanal(covmat=Harman74.cor$cov,factors=4)
factor.congruence(list(mle,pa,pc,uls,wls))
#note that the order of factors and the sign of some of factors differ 

#finally, compare the unrotated factor, ml, uls, and  wls solutions
wls <- factor.wls(Harman74.cor$cov,4,rotate="none")
pa <- factor.pa(Harman74.cor$cov,4,rotate="none")
mle <-  factanal(factors=4,covmat=Harman74.cor$cov,rotation="none")
uls <- factor.minres(Harman74.cor$cov,4,rotate="none")

factor.congruence(list(mle,pa,uls,wls))
#note that the order of factors and the sign of some of factors differ 



#an example of where the ML and PA and MR models differ is found in Thurstone.33.
#compare the first two factors with the 3 factor solution 
Thurstone.33 <- as.matrix(Thurstone.33)
mle2 <- factanal(covmat=Thurstone.33,factors=2,rotation="none")
mle3 <- factanal(covmat=Thurstone.33,factors=3 ,rotation="none")
pa2 <- factor.pa(Thurstone.33,2,rotate="none")
pa3 <- factor.pa(Thurstone.33,3,rotate="none")
mr2 <- fa(Thurstone.33,2,rotate="none")
mr3 <- fa(Thurstone.33,3,rotate="none")
factor.congruence(list(mle2,mle3,pa2,pa3,mr2,mr3))
}

\keyword{ multivariate }
\keyword{ models}