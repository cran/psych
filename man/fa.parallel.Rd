\name{fa.parallel}
\alias{fa.parallel}
\alias{fa.parallel.poly}
\alias{plot.poly.parallel}
\title{Scree plots of data or correlation matrix compared to random ``parallel" matrices }
\description{One way to determine the number of factors or components in a data matrix or a correlation matrix is to examine the ``scree" plot of the successive eigenvalues.  Sharp breaks in the plot suggest the appropriate number of components or factors to extract.  ``Parallel" analyis is an alternative technique that compares the scree of factors of the observed data with that of a random data matrix of the same size as the original. fa.parallel.poly does this for tetrachoric or polychoric analyses.
}
\usage{
fa.parallel(x, n.obs = NULL,fm="minres", fa="both", 
       main = "Parallel Analysis Scree Plots",
        n.iter=20,error.bars=FALSE,SMC=FALSE,ylabel=NULL,show.legend=TRUE)
fa.parallel.poly(x ,n.iter=10,SMC=TRUE,  fm = "minres",correct=TRUE,sim=FALSE,
       fa="both",global=TRUE) 
\method{plot}{poly.parallel}(x,show.legend=TRUE,fa="both",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or data matrix of scores.  If the matrix is square, it is assumed to be a correlation matrix.  Otherwise, correlations (with pairwise deletion) will be found }
  \item{n.obs}{n.obs=0 implies a data matrix/data.frame.  Otherwise, how many cases were used to find the correlations. }
  \item{fm}{What factor method to use. (minres, ml, uls, wls, gls, pa) See  \code{\link{fa}} for details.}
 \item{fa}{show the eigen values for a principal components (fa="pc") or a principal axis factor analysis (fa="fa") or both principal components and principal factors (fa="both")}
 \item{main}{ a title for the analysis }
 \item{n.iter}{Number of simulated analyses to perform}
 \item{correct}{For tetrachoric correlations, should a correction for continuity be applied. (See tetrachoric)}
 \item{sim}{In the case of polychoric or tetrachoric data, in addition to randomizing the real data, should we compare the solution to random simulated data.  This will double the processing time, but basically show the same result.}
 \item{error.bars}{Should error.bars be plotted (default = FALSE)}
 \item{SMC}{SMC=TRUE finds eigen values after estimating communalities by using SMCs.  smc = FALSE finds eigen values after estimating communalities with the first factor.}
 \item{ylabel}{Label for the y axis -- defaults to ``eigen values of factors and components", can be made empty to show many graphs}
 \item{show.legend}{the default is to have a legend.  For multiple panel graphs, it is better to not show the legend}
 \item{global}{If doing polychoric analyses (fa.parallel.poly) and the number of alternatives differ across items, it is necessary to turn off the global option}
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
 
 A sad observation about parallel analysis is that it is sensitive to sample size.  That is, for large data sets, the eigen values of random data are very close to 1.  This will lead to different estimates of the number of factors as a function of sample size.  Consider factor structure of the bfi data set (the first 25 items are meant to represent a five factor model).  For samples of 200 or less, parallel analysis suggests 5 factors, but for 1000 or more, six factors and components are indicated.  This is not due to an instability of the eigen values of the real data, but rather the closer approximation to 1 of the random data as n increases.
 
 When simulating dichotomous data in fa.parallel.poly, the simulated data have the same difficulties as the original data.  This functionally means that the simulated and the resampled results will be very similar.
 
 As with many psych functions, fa.parallel has been changed to allow for multicore processing.  For running a large number of iterations, it is obviously faster to increase the number of cores to the maximum possible (using the options("mc.cores=n) command where n is determined from detectCores(). 

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

#test.data <- Harman74.cor$cov   #The 24 variable Holzinger - Harman problem
#fa.parallel(test.data,n.obs=145)
fa.parallel(Thurstone,n.obs=213)   #the 9 variable Thurstone problem

#set.seed(123)
#minor <- sim.minor(24,4,400) #4 large and 12 minor factors
#fa.parallel(minor$observed) #shows 5 factors and 4 components -- compare with 
#fa.parallel(minor$observed,SMC=FALSE) #which shows 6  and 4 components factors


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

