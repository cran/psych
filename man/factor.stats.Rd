\name{factor.stats}
\Rdversion{1.1}
\alias{factor.stats}

\title{Find various goodness of fit statistics for factor analysis and principal components }

\description{Chi square and other goodness of fit statistics are found based upon the fit of a factor or components model to a correlation matrix.  Although these statistics are normally associated with a maximum likelihood solution, they can be found for minimal residual (OLS), principal axis, or principal component solutions as well.  Primarily called from within these functions, factor.stats can be used by itself. Measures of factorial adequacy and validity follow the paper by Grice, 2001.
}
\usage{
factor.stats(r=NULL,f,phi=NULL,n.obs=NA,np.obs=NULL,alpha=.1) 
}

\arguments{
  \item{r}{A correlation matrix or a data frame of raw data}
  \item{f}{A factor analysis loadings matrix or the output from a factor or principal components analysis.  In which case the r matrix need not be specified.}
  \item{phi}{A factor intercorrelation matrix if the factor solution was oblique.}
  \item{n.obs}{The number of observations for the correlation matrix.  If not specified, and a correlation matrix is used, chi square will not be reported. Not needed if the input is a data matrix.}
  \item{np.obs}{The pairwise number of subjects for each pair in the correlation matrix.  This is used for finding observed chi square.}
  \item{alpha}{alpha level of confidence intervals for RMSEA}
}
\details{Combines the goodness of fit tests used in factor.pa, factor.minres, and principal into one function.  If the matrix is singular, will smooth the correlation matrix before finding the fit functions. Now will find the RMSEA (root mean square error of approximation) and the alpha confidence intervals similar to a SEM function.  Also reports the root mean square residual.

Chi square is found two ways.  The first (STATISTIC) applies the goodness of fit test from Maximum Likelihood objective function (see below).  This assumes multivariate normality.  The second is the empirical chi square based upon the observed residual correlation matrix and the observed sample size for each correlation.  This is found by summing the squared residual correlations time the sample size.  
}
\value{ 
\item{fit}{How well does the factor model reproduce the correlation matrix. (See \code{\link{VSS}}, \code{\link{ICLUST}}, and \code{\link{principal}} for this fit statistic.}
  \item{fit.off}{how well are the off diagonal elements reproduced?  This is just 1 - the  relative magnitude of the squared off diagonal residuals to the squared off diagonal original values.}
  \item{dof}{Degrees of Freedom for this model. This is the number of observed correlations minus the number of independent parameters.  Let n=Number of items, nf = number of factors then
  \cr
  \eqn{dof = n * (n-1)/2 - n * nf + nf*(nf-1)/2}{dof = n * (n-1)/2 - n * nf + nf*(nf-1)/2}}
  \item{objective}{value of the function that is minimized by maximum likelihood procedures.  This is reported for comparison purposes and as a way to estimate chi square goodness of fit.  The objective function is 
  \cr
  \eqn{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}. }
  
  \item{STATISTIC}{If the number of observations is specified or found, this is a chi square based upon the objective function, f.  Using the formula from \code{\link{factanal}}(which seems to be Bartlett's test) :
 \cr
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{ chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } 
 
 Note that this is different from the chi square reported by the sem package which seems to use 
 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 )* f } 
 
 }
 
 \item{PVAL}{If n.obs > 0, then what is the probability of observing a chisquare this large or larger?}
 
 \item{Phi}{If oblique rotations (using oblimin from the GPArotation package or promax) are requested, what is the interfactor correlation.}

\item{R2}{The multiple R square between the factors and factor score estimates, if they were to be found. (From Grice, 2001)}

\item{r.scores}{The correlations of the factor score estimates, if they were to be found.}

\item{weights}{The beta weights to find the factor score estimates}

\item{valid}{The validity coffiecient of course coded (unit weighted) factor score estimates (From Grice, 2001)}

\item{score.cor}{The correlation matrix of course coded (unit weighted) factor score estimates, if they were to be found, based upon the loadings matrix.}

\item{RMSEA}{The Root Mean Square Error of Approximation and the alpha confidence intervals. Based upon the chi square non-centrality parameter.
This is found as \eqn{\sqrt{f/dof - 1(/-1)}}{sqrt((f/df - 1(N-1))} }
\item{rms}{The empirically found square root of the squared residuals.  This does not require sample size to be specified nor does it make assumptions about normality.}
\item{crms}{While the rms uses the number of correlations to find the average, the crms uses the number of degrees of freedom.  Thus, there is a penalty for having too complex a model.}
}
\author{William Revelle}


\references{
Grice, James W.,2001, Computing and evaluating factor scores,  Psychological Methods, 6,4, 430-450.

 } 
 
\seealso{
\code{\link{fa}} with fm="pa" for principal axis factor analysis, \code{\link{fa}} with fm="minres" for minimum residual factor analysis (default).  \code{\link{factor.pa}} also does principal axis factor analysis, but is deprecated, as is \code{\link{factor.minres}} for minimum residual factor analysis. See \code{\link{principal}} for principal components.
}
\examples{
v9 <- sim.hierarchical()
f3 <- fa(v9,3)
factor.stats(v9,f3,n.obs=500)
f3o <- fa(v9,3,fm="pa",rotate="Promax")
factor.stats(v9,f3o,n.obs=500)


}

\keyword{ multivariate }
\keyword{ models}% __ONLY ONE__ keyword per line
