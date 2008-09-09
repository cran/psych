\name{smc}
\alias{smc}
\alias{SMC}

\title{Find the Squared Multiple Correlation (SMC) of each variable with the remaining variables in a matrix}
\description{The squared multiple correlation of a variable with the remaining variables in a matrix is sometimes used as initial estimates of the communality of a variable.}

\usage{
smc(R)
}

\arguments{
  \item{R}{ A correlation matrix or a dataframe.  In the latter case, correlations are found.}
}

\value{a vector of squared multiple correlations. 
}

\author{ William Revelle }

\seealso{  \code{\link{mat.regress}}, \code{\link{factor.pa}} }
\examples{
R <- make.hierarchical()
round(smc(R),2)
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

