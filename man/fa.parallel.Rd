\name{fa.parallel}
\alias{fa.parallel}
\alias{fa.parallel.poly}
\alias{plot.poly.parallel}
\title{Scree plots of data or correlation matrix compared to random ``parallel" matrices }
\description{One way to determine the number of factors or components in a data matrix or a correlation matrix is to examine the ``scree" plot of the successive eigenvalues.  Sharp breaks in the plot suggest the appropriate number of components or factors to extract.  ``Parallel" analyis is an alternative technique that compares the scree of factors of the observed data with that of a random data matrix of the same size as the original. fa.parallel.poly does this for tetrachoric or polychoric analyses.
}
\usage{
fa.parallel(x, n.obs = NULL,fm="minres", fa="both", main = "Parallel Analysis Scree Plots",n.iter=20,error.bars=FALSE,SMC=FALSE,ylabel=NULL,show.legend=TRUE)
fa.parallel.poly(x ,n.iter=10,SMC=TRUE,  fm = "minres") 
plot.poly.parallel(x,show.legend=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or data matrix of scores.  If the matrix is square, it is assumed to be a correlation matrix.  Otherwise, correlations (with pairwise deletion) will be found }
  \item{n.obs}{n.obs=0 implies a data matrix/data.frame.  Otherwise, how many cases were used to find the correlations. }
  \item{fm}{What factor method to use. (minres, ml, uls, wls, gls, pa) See  \code{\link{fa}} for details.}
 \item{fa}{show the eigen values for a principal components (fa="pc") or a principal axis factor analysis (fa="fa") or both principal components and principal factors (fa="both")}
 \item{main}{ a title for the analysis }
 \item{n.iter}{Number of simulated analyses to perform}
 \item{error.bars}{Should error.bars be plotted (default = FALSE)}
 \item{SMC}{SMC=TRUE finds eigen values after estimating communalities by using SMCs.  smc = FALSE finds eigen values after estimating communalities with the first factor.}
 \item{ylabel}{Label for the y axis -- defaults to ``eigen values of factors and components", can be made empty to show many graphs}
 \item{show.legend}{the default is to have a legend.  For multiple panel graphs, it is better to not show the legend}
 \item{...}{additional plotting parameters, for plot.poly.parallel}
}
\details{Cattell's ``scree" test is one of most simple tests for the number of factors problem.  Horn's (1965) ``parallel" analysis is an equally compelling procedure.  Other procedures for determining the most optimal number of factors include finding the Very Simple Structure (VSS) criterion (\code{\link{VSS}}) and Velicer's \code{\link{MAP}} procedure (included in \code{\link{VSS}}).  fa.parallel plots the eigen values for a principal components and the factor solution (minres by default) and does the same for random matrices of the same size as the original data matrix.  For raw data, the random matrices are 1) a matrix of univariate normal data and 2) random samples (randomized across rows) of the original data.

\code{\link{fa.parallel.poly}} will do parallel analysis for polychoric and tetrachoric factors.  If the data are dichotomous, \code{\link{fa.parallel.poly}} will find tetrachoric correlations for the real and simulated data, otherwise, if the number of categories is less than 10, it will find polychoric correlations.  Note that fa.parallel.poly is much slower than fa.parallel because of the complexity of calculating the tetrachoric/polychoric correlations.

The means of (ntrials) random solutions are shown.  Error bars are usually very small and are suppressed by default but can be shown if requested.

Alternative ways to estimate the number of factors problem are discussed in the Very Simple Structure  (Revelle and Rocklin, 1979) documentation (\code{\link{VSS}}) and include Wayne Velicer's \code{\link{MAP}} algorithm (Veicer, 1976).  

Parallel analysis for factors is actually harder than it seems, for the question is what are the appropriate communalities to use.  If communalities are estimated by the Squared Multiple Correlation (SMC) \code{\link{smc}}, then the eigen values of the original data will reflect major as well as minor factors (see \code{\link{sim.minor}} to simulate such data).  Random data will not, of course, have any structure and thus the number of factors will tend to be biased upwards by the presence of the minor factors.  

By default, fa.parallel estimates the communalities based upon a one factor minres solution.  Although this will underestimate the communalities, it does seem to lead to better solutions on simulated or real (e.g., the \code{\link{bfi}} or Harman74) data sets.  

For comparability with other algorithms (e.g, the paran function in the paran package), setting smc=TRUE will use smcs as estimates of communalities. This will tend towards identifying more factors than the default option.

 Printing the results will show the eigen values of the original data that are greater than simulated values.

}
\value{
  A plot of the eigen values for the original data, ntrials of resampling of the  original data, and of a equivalent size matrix of random normal deviates.  If the data are a correlation matrix, specify the number of observations.  
  
  Also returned (invisibly) are:
  
  \item{fa.values}{The eigen values of the factor model for the real data.}
  \item{fa.sim}{The descriptive statistics of the simulated factor models.}
  \item{pc.values}{The eigen values of a principal components of the real data.}
  \item{pc.sim}{The descriptive statistics of the simulated principal components analysis.}
  \item{nfact}{Number of factors with eigen values > eigen values of random data}
  \item{ncomp}{Number of components with eigen values > eigen values of random data}
  
 
  
}
\references{


Floyd, Frank J.  and  Widaman, Keith. F (1995)  Factor analysis in the development and refinement of clinical assessment instruments. Psychological Assessment, 7(3):286-299, 1995.

Horn, John (1965) A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179-185.

Humphreys, Lloyd G. and Montanelli, Richard G. (1975), An investigation of the parallel analysis criterion for determining the number of common factors.
Multivariate Behavioral Research, 10, 193-205.


Revelle, William and Rocklin, Tom (1979) Very simple structure - alternative procedure for estimating the optimal number of interpretable factors. Multivariate Behavioral Research, 14(4):403-414.

Velicer, Wayne. (1976) Determining the number of components from the matrix of partial correlations. Psychometrika, 41(3):321-327, 1976.

 }
\author{ William Revelle }


\seealso{  \code{\link{fa}}, \code{\link{VSS}}, \code{\link{VSS.plot}},   \code{\link{VSS.parallel}}, \code{\link{sim.minor}}}
\examples{

test.data <- Harman74.cor$cov 
fa.parallel(test.data,n.obs=145)
set.seed(123)
minor <- sim.minor(24,4,400) #4 large and 12 minor factors
fa.parallel(minor$observed) #shows 4 factors  -- compare with 
fa.parallel(minor$observed,SMC=TRUE) #which shows 8 factors


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

