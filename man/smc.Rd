\name{smc}
\alias{smc}


\title{Find the Squared Multiple Correlation (SMC) of each variable with the remaining variables in a matrix}
\description{The squared multiple correlation of a variable with the remaining variables in a matrix is sometimes used as initial estimates of the communality of a variable.

SMCs are also used when estimating reliability using Guttman's lambda 6  \code{\link{guttman}} coefficient. 

The SMC is just  1 - 1/diag(R.inv) where R.inv is the inverse of R.
}
\usage{
smc(R,covar=FALSE)
}

\arguments{
  \item{R}{ A correlation matrix or a dataframe.  In the latter case, correlations are found.}
  \item{covar}{if covar = TRUE and  R is either a covariance matrix or data frame, then return the smc * variance for each item}
}

\value{a vector of squared multiple correlations.  Or, if covar=TRUE, a vector of squared multiple correlations * the item variances

If the matrix is not invertible, then a vector of 1s is returned
}

\author{ William Revelle }

\seealso{  \code{\link{mat.regress}}, \code{\link{fa}}  }
\examples{
R <- make.hierarchical()
round(smc(R),2)
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

