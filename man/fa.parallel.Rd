\name{fa.parallel}
\alias{fa.parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scree plots of data or correlation matrix compared to random ``parallel" matrices }
\description{One way to determine the number of factors or components in a data matrix or a correlation matrix is to examine the ``scree" plot of the successive eigenvalues.  Sharp breaks in the plot suggest the appropriate number of components or factors to extract.  ``Parallel" analyis is an alternative technique that compares the scree of the observed data with that of a random data matrix of the same size as the original.  
}
\usage{
fa.parallel(x, n.obs = NULL, fa="both", main = "Parallel Analysis Scree Plots")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or data matrix of scores.  If the matrix is square, it is assumed to be a correlation matrix.  Otherwise, correlations (with pairwise deletion) will be found }
  \item{n.obs}{n.obs=0 implies a data matrix/data.frame.  Otherwise, how many cases were used to find the correlations. }
 \item{fa}{show the eigen values for a principal components (fa="pc") or a principal axis factor analysis (fa="fa") or both principal components and principal factors (fa="both")}
 \item{main}{ a title for the analysis }
}
\details{Cattell's ``scree" test is one of most simple tests for the number of factors problem.  Humphreys and Montanelli's ``parallel" analysis is an equally compelling procedure.  Other procedures for determining the most optimal number of factors include finding the Very Simple Structure (VSS) criterion (\code{\link{VSS}}).  fa.parallel plots the eigen values for a principal components and principal factor solution and does the same for random matrices of the same size as the original data matrix.  For raw data, the random matrices are 1) a matrix of univariate normal data and 2) random samples (randomized across rows) of the original data.
}
\value{
  A plot of the eigen values for the original data, a resampling of the  original data, and of a equivalent size matrix of random normal deviates.  If the data are a correlation matrix, specify the number of observations.  
}
\references{ }
\author{ William Revelle }


\seealso{ \code{\link{VSS}},\code{\link{VSS.plot}}, \code{\link{VSS.parallel}} }
\examples{

test.data <- Harman74.cor$cov 
fa.parallel(test.data,n.obs=200)

fa.parallel(attitude) 
#

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

