\name{fa.poly}
\alias{factor.pa} 
\alias{factor.minres} 
\alias{factor.wls}  
\alias{fa.poly}


\title{Deprecated Exploratory Factor analysis functions.  Please use fa}
\description{After 4 years, it is time to stop using these deprecated functions!  Please see \code{\link{fa}} which includes all of the functionality of these older functions.  
}

\usage{

fa.poly(x,nfactors=1,n.obs = NA, n.iter=1, rotate="oblimin", SMC=TRUE,  missing=FALSE,
 impute="median", min.err = .001, max.iter=50, symmetric=TRUE, warnings=TRUE,
 fm="minres",alpha=.1, p =.05,scores="regression", oblique.scores=TRUE,
        weight=NULL,global=TRUE,...)  #deprecated
        
factor.minres(r, nfactors=1, residuals = FALSE, rotate = "varimax",n.obs = NA,
scores = FALSE,SMC=TRUE, missing=FALSE,impute="median",min.err = 0.001, digits = 2,
 max.iter = 50,symmetric=TRUE,warnings=TRUE,fm="minres")    #deprecated

factor.wls(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NA,
scores=FALSE,SMC=TRUE,missing=FALSE,impute="median", min.err = .001,
 digits=2,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="wls")  #deprecated
}

\arguments{
  \item{r}{deprecated.}
  \item{x}{deprecated}
  \item{nfactors}{ deprecated }
  \item{n.obs}{deprecated }
   \item{rotate}{deprecated}
 \item{n.iter}{deprecated}
  \item{residuals}{deprecated }
  \item{scores}{deprecated}
  \item{SMC}{deprecated }

  \item{missing}{deprecated}
  \item{impute}{deprecated}
  \item{max.iter}{deprecated}
  \item{symmetric}{deprecated}
  \item{warnings}{deprecated}
  \item{fm}{deprecated}
  \item{alpha}{deprecated}
  \item{p}{deprecated}
  \item{oblique.scores}{deprecated}
  \item{weight}{deprecated}
  
  \item{global}{deprecated}
  \item{digits}{deprecated}
  \item{min.err}{deprecated}
  \item{...}{deprecated}
  
 } 
 
\details{Please see the writeup for \code{\link{fa}} for all of the functionality in these older functions.
}
\value{
  \item{values }{Eigen values of the common factor solution}
  \item{e.values}{Eigen values of the original matrix}
  \item{communality}{Communality estimates for each item.  These are merely the sum of squared factor loadings for that item.}
  \item{communalities}{If using minrank factor analysis, these are the communalities reflecting the total amount of common variance.  They will exceed the communality (above) which is the model estimated common variance. }
  \item{rotation}{which rotation was requested?}
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
 
 \item{Phi}{If oblique rotations (using oblimin from the GPArotation package or promax) are requested, what is the interfactor correlation.}
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
}

\keyword{ multivariate }
\keyword{ models}