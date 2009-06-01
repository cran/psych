\name{factor.stats}
\Rdversion{1.1}
\alias{factor.stats}
\alias{factor.scores}

\title{Find various goodness of fit statistics for factor analysis and principal components }

\description{Chi square and other goodness of fit statistics are found based upon the fit of a factor or components model to a correlation matrix.  Although these statistics are normally associated with a maximum likelihood solution, they can be found for minimal residual (OLS), principal axis, or principal component solutions as well.  Primarily called from within these functions, factor.stats can be used by itself. Measures of factorial adequacy and validity follow the paper by Grice, 2001.
}
\usage{
factor.stats(r, f, phi=NULL, n.obs = NA)
factor.scores(x,f)
}

\arguments{
  \item{r}{A correlation matrix or a data frame of raw data}
  \item{x}{A data frame or matrix of raw data}
  \item{f}{A factor analysis loadings matrix or the output from a factor or principal components analysis}
  \item{phi}{A factor intercorrelation matrix if the factor solution was oblique}
  \item{n.obs}{The number of observations for the correlation matrix.  If not specified, and a correlation matrix is used, chi square will not be reported. Not needed if the input is a data matrix.
}
}
\details{Combines the goodness of fit tests used in factor.pa, factor.minres, and principal into one function.  If the matrix is singular, will smooth the correlation matrix before finding the fit functions.
}
\value{ 
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

\item{R2}{The multiple R square between the factors and factor score estimates, if they were to be found. (From Grice, 2001)}

\item{r.scores}{The correlations of the factor score estimates, if they were to be found.}

\item{weights}{The beta weights to find the factor score estimates}

\item{valid}{The validity coffiecient of course coded (unit weighted) factor score estimates (From Grice, 2001)}

\item{score.cor}{The correlation matrix of course coded (unit weighted) factor score estimates, if they were to be found, based upon the loadings matrix.}

}
\author{William Revelle}


\references{
Grice, James W.,2001, Computing and evaluating factor scores,  Psychological Methods, 6,4, 430-450.

 } 
 
\seealso{
\code{\link{factor.pa}} for principal axis factor analysis, \code{\link{factor.minres}} for minimum residual factor analysis, and \code{\link{principal}} for principal components.
}
\examples{
v9 <- sim.hierarchical()
f3 <- factor.minres(v9,3)
factor.stats(v9,f3,n.obs=500)
f3o <- factor.pa(v9,3,rotate="Promax")
factor.stats(v9,f3o,n.obs=500)


}

\keyword{ multivariate }
\keyword{ models}% __ONLY ONE__ keyword per line
