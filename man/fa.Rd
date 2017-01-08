\name{fa}
\alias{fa}
\alias{fac}



\title{Exploratory Factor analysis using MinRes (minimum residual) as well as EFA by Principal Axis, Weighted Least Squares or Maximum Likelihood }
\description{Among the many ways to do latent variable exploratory factor analysis (EFA), one of the better is to use Ordinary Least Squares (OLS) to find the minimum residual (minres) solution. This produces solutions very similar to maximum likelihood even for badly behaved matrices. A variation on minres is to do weighted least squares (WLS). Perhaps the most conventional technique is principal axes (PAF).  An eigen value decomposition of a correlation matrix is done and then the communalities for each variable are estimated by the first n factors. These communalities are entered onto the diagonal and the procedure is repeated until the sum(diag(r)) does not vary.   Yet another estimate procedure is maximum likelihood. For well behaved matrices, maximum likelihood factor analysis (either in the fa or in the factanal function) is probably preferred.  Bootstrapped confidence intervals of the loadings and interfactor correlations are found by fa  with n.iter > 1. }

\usage{
fa(r,nfactors=1,n.obs = NA,n.iter=1, rotate="oblimin", scores="regression", 
residuals=FALSE, SMC=TRUE, covar=FALSE,missing=FALSE,impute="median",
min.err = 0.001,  max.iter = 50,symmetric=TRUE, warnings=TRUE, fm="minres",
 alpha=.1,p=.05,oblique.scores=FALSE,np.obs,use="pairwise",cor="cor",weight=NULL,...)

fac(r,nfactors=1,n.obs = NA, rotate="oblimin", scores="tenBerge", residuals=FALSE,
 SMC=TRUE, covar=FALSE,missing=FALSE,impute="median",min.err = 0.001, 
max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1,
oblique.scores=FALSE,np.obs,use="pairwise",cor="cor",weight=NULL,...)

}

\arguments{
  \item{r}{A correlation or covariance matrix or a raw data matrix. If raw data, the correlation matrix will be found using pairwise deletion. If covariances are supplied, they will be converted to correlations unless the covar option is TRUE.}
  \item{nfactors}{ Number of factors to extract, default is 1 }
  \item{n.obs}{Number of observations used to find the correlation matrix if using a correlation matrix.  Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.}
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
  \item{fm}{factoring method  fm="minres" will do a minimum residual (OLS), fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample size weighted chi square when treating pairwise correlations with different number of subjects per pair. fm ="minrank" will do a minimum rank factor analysis.}
  \item{alpha}{alpha level for the confidence intervals for RMSEA}
  \item{p}{if doing iterations to find confidence intervals, what probability values should be found for the confidence intervals}
  \item{oblique.scores}{When factor scores are found, should they be based on the structure matrix (default) or the pattern matrix (oblique.scores=TRUE).  } 
  \item{weight}{If not NULL, a vector of length n.obs that contains weights for each observation. The NULL case is equivalent to all cases being weighted 1.}
  \item{use}{How to treat missing data, use="pairwise" is the default".  See cor for other options.}
  \item{cor}{How to find the correlations: "cor" is Pearson", "cov" is covariance, 
  "tet" is tetrachoric, "poly" is polychoric, "mixed" uses mixed cor for a mixture of tetrachorics, polychorics, Pearsons, biserials, and polyserials, Yuleb is Yulebonett, Yuleq and YuleY are the obvious Yule coefficients as appropriate}

  \item{...}{additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax}
  
 } 
 
\details{Factor analysis is an attempt to approximate a correlation or covariance matrix with one of lesser rank.  The basic model is that \eqn{_nR_n \approx _{n}F_{kk}F_n'+ U^2}{nRn = nFk kFn' + U2} where k is much less than n. There are many ways to do factor analysis, and maximum likelihood procedures are probably the most preferred (see \code{\link{factanal}} ).  The existence of uniquenesses is what distinguishes factor analysis from principal components analysis (e.g., \code{\link{principal}}). If variables are thought to represent a ``true" or latent part then factor analysis provides an estimate of the correlations with the latent factor(s) representing the data.  If variables are thought to be measured without error, then principal components provides the most parsimonious description of the data.  

The fa function will do factor analyses using one of five different algorithms: minimum residual (minres), principal axes, weighted least squares, minimum rank, or maximum likelihood.

Principal axes factor analysis has a long history in exploratory analysis and is a straightforward procedure.  Successive eigen value decompositions are done on a correlation matrix with the diagonal replaced with  diag (FF') until \eqn{\sum(diag(FF'))} does not change (very much).  The current limit of max.iter =50 seems to work for most problems, but the Holzinger-Harmon 24 variable problem needs about 203 iterations to converge for a 5 factor solution.  

Not all factor programs that do principal axes do iterative solutions.  The example from the SAS manual (Chapter 26) is such a case. To achieve that solution, it is necessary to specify that the max.iterations = 1.  Comparing that solution to an iterated one (the default) shows that iterations improve the solution.  In addition, fm="minres" or fm="mle" produces even better solutions for this example.  

Principal axes may be used in cases when maximum likelihood solutions fail to converge, although fm="minres" will also do that and tends to produce better (smaller residuals) solutions.

The fm="minchi" option is a variation on the "minres" (ols) solution and minimizes the sample size weighted residuals rather than just the residuals. This was developed to handle the problem of data that Massively Missing Completely at Random (MMCAR) which a condition that happens in the SAPA project.


A problem in factor analysis is to find the best estimate of the original communalities.  Using the Squared Multiple Correlation (SMC) for each variable will underestimate the communalities, using 1s will over estimate.  By default, the SMC estimate is used.  In either case, iterative techniques will tend to converge on a stable solution. If, however, a solution fails to be achieved, it is useful to try again using ones (SMC =FALSE).  Alternatively, a vector of starting values for the communalities may be specified by the SMC option.

The iterated principal axes algorithm does not attempt to find the best (as defined by a maximum likelihood criterion) solution, but rather one that converges rapidly using successive eigen value decompositions.  The maximum likelihood criterion of fit and the associated chi square value are reported, and will be worse than that found using maximum likelihood procedures.

The minimum residual (minres) solution is an unweighted least squares solution that takes a slightly different approach.  It uses the \code{\link{optim}} function and adjusts the diagonal elements of the correlation matrix to mimimize the squared residual when the factor model is the eigen value decomposition of the reduced matrix.  MINRES and PA will both work when ML will not, for they can be used when the matrix is singular. At least on a number of test cases, the MINRES solution is slightly more similar to the ML solution than is the PA solution.  To a great extent, the minres and wls solutions follow ideas in the \code{\link{factanal}} function.  

The weighted least squares (wls) solution weights the residual matrix by 1/ diagonal of the inverse of the correlation matrix.  This has the effect of weighting items with low communalities more than those with high communalities. 

The generalized least squares (gls) solution weights the residual matrix by the inverse of the correlation matrix.  This has the effect of weighting those variables with low communalities even more than those with high communalities.

The maximum likelihood solution takes yet another approach and finds those communality values that minimize the chi square goodness of fit test.  The fm="ml" option provides a maximum likelihood solution following the procedures used in \code{\link{factanal}} but does not provide all the extra features of that function.

The minimum rank factor model (MRFA) roughly follows ideas by Shapiro and Ten Berge (2002) and Ten Berge and Kiers (1991).  It makes use of the glb.algebraic procedure contributed by Andreas Moltner.  MRFA attempts to extract factors such that the residual matrix is still positive semi-definite.  This version is still being tested and feedback is most welcome.

Test cases comparing the output to SPSS suggest that the PA algorithm matches what SPSS calls uls, and that the wls solutions are equivalent in their fits. The wls and gls solutions have slightly larger eigen values, but slightly worse fits of the off diagonal residuals than do the minres or maximum likelihood solutions.  Comparing the results to the examples in Harman (76), the PA solution with no iterations matches what Harman calls Principal Axes (as does SAS), while the iterated PA solution matches his minres solution.  The minres solution found in psych tends to have slightly smaller off diagonal residuals (as it should) than does the iterated PA solution.   

Although for items, it is typical to find factor scores by scoring the salient items (using, e.g., \code{\link{scoreItems}}) factor scores can be estimated by regression as well as several other means. There are multiple approaches that are possible (see Grice, 2001) and the one taken here was developed by tenBerge et al.  (see \code{\link{factor.scores}}. The alternative, which will match factanal is to find the scores using regression --  Thurstone's least squares regression where the weights are found by   \eqn{W = R^(-1)S}{W = inverse(R)S} where R is the correlation matrix of the variables ans S is the structure matrix.  Then, factor scores are just \eqn{Fs = X W}{Fs = X W}.

In the oblique case, the factor loadings are referred to as Pattern coefficients and are related to the Structure coefficients by \eqn{S = P \Phi}{S = P Phi} and thus \eqn{P = S \Phi^{-1}}{P = S Phi^-1}.  When estimating factor scores, \code{\link{fa}} and  \code{\link{factanal}} differ in that \code{\link{fa}} finds the factors from the Structure matrix while \code{\link{factanal}} seems to do it from the Pattern matrix.  Thus, although in the orthogonal case, fa and factanal agree perfectly in their factor score estimates, they do not agree in the case of oblique factors.  Setting oblique.scores = TRUE  will produce factor score estimate that match those of \code{\link{factanal}}.

It is sometimes useful to extend the factor solution to variables that were not factored.  This may be done using \code{\link{fa.extension}}. Factor extension is typically done in the case where some variables were not appropriate to factor, but factor loadings on the original factors are still desired.

For dichotomous items or polytomous items, it is recommended to analyze the \code{\link{tetrachoric}} or \code{\link{polychoric}} correlations rather than the Pearson correlations. This may be done by specifying cor="poly" or cor="tet" or cor="mixed" if the data have a mixture of dichotomous, polytomous, and continous variables.  

Analysis of dichotomous or polytomous data may also be done by using  \code{\link{irt.fa}} or simply setting the cor="poly" option.  In the first case,  the factor analysis results are reported in Item Response Theory (IRT) terms, although the original factor solution is returned in the results. In the later case, a typical factor loadings matrix is returned, but the tetrachoric/polychoric correlation matrix and item statistics are saved for reanalysis by \code{\link{irt.fa}}. (See also the \code{\link{mixed.cor}} function to find correlations from a mixture of continuous, dichotomous, and polytomous items.)


Of the various rotation/transformation options, varimax, Varimax, quartimax, bentlerT, geominT, and bifactor do orthogonal rotations. Promax  transforms obliquely with a target matix equal to the varimax solution. oblimin, quartimin, simplimax,  bentlerQ,  geominQ and biquartimin are oblique transformations. Most of these are just calls to  the GPArotation package. The ``cluster'' option does a targeted rotation to a structure defined by the cluster representation of a varimax solution.  With the optional "keys" parameter, the "target" option will rotate to a target supplied as a keys matrix. (See \code{\link{target.rot}}.)

Two additional target rotation options are available through calls to GPArotation.  These are the targetQ (oblique) and targetT (orthogonal) target rotations of Michael Browne.  See \code{\link{target.rot}} for more documentation. 

The "bifactor" rotation implements the Jennrich and Bentler (2011) bifactor rotation by calling the GPForth function in the GPArotation package and using two functions adapted from  the MatLab code of Jennrich and Bentler.

There are two varimax rotation functions.  One, Varimax, in the GPArotation package does not by default apply Kaiser normalization.  The other, varimax, in the stats package, does.  It appears that the two rotation functions produce slightly different results even when normalization is set. For consistency with the other rotation functions, Varimax is probably preferred.

The rotation matrix  (rot.mat) is returned from all of these options. This is the inverse of the Th (theta?) object returned by the GPArotation package.  The correlations of the factors may be found by \eqn{\Phi = \theta' \theta}{Phi = Th' Th}

There are two ways to handle dichotomous or polytomous responses: \code{\link{fa}} with the cor="poly" option which will return the tetrachoric or polychoric correlation matrix, as well as the normal factor analysis output, and \code{\link{irt.fa}} which returns a two parameter irt analysis as well as the normal fa output. 

When factor analyzing items with dichotomous or polytomous responses, the \code{\link{irt.fa}} function provides an Item Response Theory representation of the factor output. The factor analysis results are available, however, as an object in the irt.fa output. 

\code{\link{fa.poly}} is deprecated, for its functioning is matched by setting cor="poly".  It will produce normal factor analysis output but also will save the polychoric matrix  (rho) and items difficulties (tau) for subsequent irt analyses.  \code{\link{fa.poly}} will,  by default, find factor scores if the data are available.  The correlations are found using either \code{\link{tetrachoric}} or \code{\link{polychoric}} and then this matrix is factored.  Weights from the factors are then applied to the original data to estimate factor scores.
  

The function \code{\link{fa}} will repeat the analysis n.iter times on a bootstrapped sample of the data (if they exist) or of a simulated data set based upon the observed correlation matrix.  The mean estimate and standard deviation of the estimate are returned and will print the original factor analysis as well as the alpha level confidence intervals for the estimated coefficients.  The bootstrapped solutions are rotated towards the original solution using target.rot. The factor loadings are z-transformed, averaged and then back transformed. This leads to an error in the case of Heywood cases.  The probably better alternative is to just find the mean bootstrapped value and find the confidence intervals based upon the observed range of values. The default is to have n.iter =1 and thus not do bootstrapping.

If using polytomous or dichotomous items, it is perhaps more useful is to find the Item Response Theory parameters equivalent to the factor loadings reported in fa.poly by using the \code{\link{irt.fa}} function.  

Some correlation matrices that arise from using pairwise deletion or from tetrachoric or polychoric matrices will not be proper.  That is, they will not be positive semi-definite (all eigen values >= 0).  The \code{\link{cor.smooth}} function will adjust correlation matrices (smooth them) by making all negative eigen values slightly greater than 0, rescaling the other eigen values to sum to the number of variables, and then recreating the correlation matrix.  See \code{\link{cor.smooth}} for an example of this problem using the \code{\link{burt}} data set.

For those who like SPSS type output, the measure of factoring adequacy known as the Kaiser-Meyer-Olkin \code{\link{KMO}} test may be found from the correlation matrix or data matrix using the \code{\link{KMO}} function.  Similarly, the Bartlett's test of Sphericity may be found using the  \code{\link{cortest.bartlett}} function.

For those who want to have an object of the variances accounted for, this is returned invisibly by the print function.  (e.g., p <- print(fa(ability))$Vaccounted )

The output from the print.psych.fa function displays the factor loadings (from the pattern matrix, the h2 (communalities) the u2 (the uniquenesses), com (the complexity of the factor loadings for that variable (see below).  In the case of an orthogonal solution, h2 is merely the row sum of the squared factor loadings. But for an oblique solution, it is the row sum of the orthogonal factor loadings (remember, that rotations or transformations do not change the communality).  

}
\value{
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

\item{weights}{The beta weights to find the factor score estimates.  These are also used by the \code{\link{predict.psych}} function to find predicted factor scores for new cases.  }
\item{scores}{The factor scores as requested.  Note that these scores reflect the choice of the way scores should be estimated (see scores in the input).  That is, simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms (see \code{\link{factor.scores}}).  The correlation between factor score estimates (r.scores) is based upon using the regression/Thurstone approach.  The actual correlation between scores will reflect the rotation algorithm chosen and may be found by correlating those scores.}

\item{valid}{The validity coffiecient of course coded (unit weighted) factor score estimates (From Grice, 2001)}

\item{score.cor}{The correlation matrix of course coded (unit weighted) factor score estimates, if they were to be found, based upon the loadings matrix rather than the weights matrix.  }

\item{rot.mat}{The rotation matrix as returned from GPArotation.}
 }



\references{Gorsuch, Richard, (1983) Factor Analysis. Lawrence Erlebaum Associates.

Grice, James W.  (2001), Computing and evaluating factor scores. Psychological Methods, 6, 430-450

Harman, Harry and Jones, Wayne (1966) Factor analysis by minimizing residuals (minres), Psychometrika, 31, 3, 351-368.

Hofmann, R. J. ( 1978 ) . Complexity and simplicity as objective indices descriptive of factor solutions. Multivariate Behavioral Research, 13, 247-250.

Pettersson E, Turkheimer E. (2010) Item selection, evaluation, and simple structure in personality data. Journal of research in personality, 44(4), 407-420.

Revelle, William. (in prep) An introduction to psychometric theory with applications in R. Springer.  Working draft available at \url{http://personality-project.org/r/book/} 

Shapiro, A. and ten Berge, Jos M. F, (2002) Statistical inference of minimum rank factor analysis.  Psychometika, (67) 79-84.

ten Berge, Jos M. F. and Kiers, Henk A. L. (1991). A numerical approach to the approximate and the exact minimum rank of a covariance matrix. Psychometrika, (56) 309-315. 
}
\author{ William Revelle }

\note{Thanks to Erich Studerus for some very helpful suggestions about various rotation and factor scoring algorithms, and to  Gumundur Arnkelsson  for suggestions about factor scores for singular matrices.

The fac function is the original fa function which is now called by fa repeatedly to get confidence intervals. 

SPSS will sometimes use a Kaiser normalization before rotating.  This will lead to different solutions than reported here.  To get the Kaiser normalized loadings, use \code{\link{kaiser}}.

The communality for a variable is the amount of variance accounted for by all of the factors.  That is to say, for orthogonal factors, it is the sum of the squared factor loadings (rowwise).  The communality is insensitive to rotation.  However, if an oblique solution is found, then the communality is not the sum of squared pattern coefficients. In both cases (oblique or orthogonal) the communality is the diagonal of the reproduced correlation matrix where \eqn{_nR_n =  _{n}P_{k  k}\Phi_{k k}P_n'}{nRn = nPk k\Phi k kPn' }  where P is the pattern matrix and \eqn{\Phi} is the factor intercorrelation matrix. This is the same, of course to multiplying the pattern by the structure: \eqn{R = P S'} {R = PS'} where the Structure matrix is \eqn{S = \Phi P}{S = Phi P}. Similarly, the eigen values are the diagonal of the product \eqn{ _k\Phi_{kk}P'_{nn}P_{k}
}{\Phi_{k k}P'_nnP_k}. 
 

A frequently asked question is why are the factor names of the rotated solution not in ascending order?  That is, for example, if factoring the 25 items of the bfi, the factor names are MR2   MR3   MR5   MR1   MR4, rather than the seemingly more logical  "MR1" "MR2" "MR3" "MR4" "MR5".  This is for pedagogical reasons, in that factors as extracted are orthogonal and are in order of amount of variance accounted for.  But when rotated (orthogonally) or transformed (obliquely) the simple structure solution does not preserve that order.  The factor names are, of course, arbitrary, and are kept with the original names to show the effect of rotation/transformation.  To give them names associated with their ordinal position, simply paste("F", 1:nf, sep="") where nf is the number of factors. See the last example.

The print function for the fa output will return (invisibly) an object (Vaccounted) that matches the printed output for the variance accounted for by each factor, as well as the cumulative variance, and the percentage of variance accounted for by each factor. 

Correction to documentation: as of September, 2014, the oblique.scores option is correctly explained. (It had been backwards.)  The default (oblique.scores=FALSE) finds scores based upon the Structure matrix, while oblique.scores=TRUE finds them based upon the pattern matrix.  The latter case matches factanal.  This error was detected by Mark Seeto. }

\seealso{ \code{\link{principal}} for principal components analysis (PCA).  PCA will give very similar solutions to factor analysis when there are many variables.  The differences become more salient as the number variables decrease.  The PCA and FA models are actually very different and should not be confused.  One is a model of the observed variables, the other is a model of latent variables.  

\code{\link{irt.fa}} for Item Response Theory analyses using factor analysis, using the two parameter IRT equivalent of loadings and difficulties.

\code{\link{VSS}} will produce the Very Simple Structure (VSS) and MAP criteria for the number of factors, \code{\link{nfactors}} to compare many different factor criteria.

\code{\link{ICLUST}} will do a hierarchical cluster analysis alternative to factor analysis or principal components analysis.

 \code{\link{predict.psych}} to find predicted scores based upon new data, \code{\link{fa.extension}} to extend the factor solution to new variables, \code{\link{omega}} for hierarchical factor analysis with one general factor. 
 code{\link{fa.multi}} for hierarchical factor analysis with an arbitrary number of higher order factors. 
 

\code{\link{fa.sort}} will sort the factor loadings into echelon form. \code{\link{fa.organize}} will reorganize the factor pattern matrix into any arbitrary order of factors and items.  

\code{\link{KMO}} and \code{\link{cortest.bartlett}} for various tests that some people like. 

\code{\link{factor2cluster}} will prepare unit weighted scoring keys of the factors that can be used with \code{\link{scoreItems}}.

\code{\link{fa.lookup}} will print the factor analysis loadings matrix along with the item ``content" taken from a dictionary of items.  This is useful when examining the meaning of the factors.  

\code{\link{anova.psych}} allows for testing the difference between two (presumably nested) factor models .




}
\examples{
#using the Harman 24 mental tests, compare a principal factor with a principal components solution
pc <- principal(Harman74.cor$cov,4,rotate="varimax")   #principal components
pa <- fa(Harman74.cor$cov,4,fm="pa" ,rotate="varimax")  #principal axis 
uls <- fa(Harman74.cor$cov,4,rotate="varimax")          #unweighted least squares is minres
wls <- fa(Harman74.cor$cov,4,fm="wls")       #weighted least squares

#to show the loadings sorted by absolute value
print(uls,sort=TRUE)

#then compare with a maximum likelihood solution using factanal
mle <- factanal(covmat=Harman74.cor$cov,factors=4)
factor.congruence(list(mle,pa,pc,uls,wls))
#note that the order of factors and the sign of some of factors may differ 

#finally, compare the unrotated factor, ml, uls, and  wls solutions
wls <- fa(Harman74.cor$cov,4,rotate="none",fm="wls")
pa <- fa(Harman74.cor$cov,4,rotate="none",fm="pa")
minres <-  factanal(factors=4,covmat=Harman74.cor$cov,rotation="none")
mle <- fa(Harman74.cor$cov,4,rotate="none",fm="mle")
uls <- fa(Harman74.cor$cov,4,rotate="none",fm="uls")
factor.congruence(list(minres,mle,pa,wls,uls))
#in particular, note the similarity of the mle and min res solutions
#note that the order of factors and the sign of some of factors may differ 



#an example of where the ML and PA and MR models differ is found in Thurstone.33.
#compare the first two factors with the 3 factor solution 
Thurstone.33 <- as.matrix(Thurstone.33)
mle2 <- fa(Thurstone.33,2,rotate="none",fm="mle")
mle3 <- fa(Thurstone.33,3 ,rotate="none",fm="mle")
pa2 <- fa(Thurstone.33,2,rotate="none",fm="pa")
pa3 <- fa(Thurstone.33,3,rotate="none",fm="pa")
mr2 <- fa(Thurstone.33,2,rotate="none")
mr3 <- fa(Thurstone.33,3,rotate="none")
factor.congruence(list(mle2,mr2,pa2,mle3,pa3,mr3))

#f5 <- fa(bfi[1:25],5)
#f5  #names are not in ascending numerical order (see note)
#colnames(f5$loadings) <- paste("F",1:5,sep="")
#f5

#Get the variance accounted for object from the print function
p <- print(mr3)
round(p$Vaccounted,2)
}

\keyword{ multivariate }
\keyword{ models}