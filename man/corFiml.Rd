\name{corFiml}
\alias{corFiml}

\title{Find a Full Information Maximum Likelihood (FIML) correlation or covariance matrix from a data matrix with missing data }
\description{Makes use of functions adapted from the lavaan package to find FIML covariance/correlation matrices.  FIML can be much slower than the normal pairwise deletion option of cor, but provides slightly more precise estimates.
}
\usage{
corFiml(x, covar = FALSE,show=FALSE) 
}
\arguments{
  \item{x}{A data.frame or data matrix}
  \item{covar}{By default, just return the correlation matrix.  If covar is TRUE, return a list containing the covariance matrix and the ML fit function.}
  \item{show}{If show=TRUE, then just show the patterns of missingness, but don't do the FIML.  Useful for understanding the process of fiml.}
}
\details{In the presence of missing data, Full Information Maximum Likelihood (FIML) is an alternative to simply using the pairwise correlations. The implementation in the lavaan package for structural equation modeling has been adapted for the simpler case of just finding the correlations or covariances.  

The pairwise solution for any pair of variables is insensitive to other variables included in the matrix.  On the other hand, the ML solution depends upon the entire set of items being correlated.  This will lead to slightly different solutions for different subsets of variables.  

The basic FIML algorithm is to find the pairwise ML solution for covariances and means for every pattern of missingness and then to weight the solution by the size of every unique pattern of missingness.  

}
\value{
\item{cor}{The correlation matrix found using FIML}
\item{cov}{The covariance matrix found using FIML}
\item{fx}{The ML fit function}
}
\author{Wiliam Revelle}
\note{
The functions used in lavaan are not exported and so have been copied (and simplified) to the psych package.

}


\seealso{ To use the resulting correlations, see \code{\link{fa}}.  To see the pairwise pattern of missingness, see \code{\link{count.pairwise}}. 
}
\examples{ 
rML <- corFiml(bfi[20:27])
rpw <- cor(bfi[20:27],use="pairwise") 
round(rML - rpw,3)
mp <- corFiml(bfi[20:27],show=TRUE)
mp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
