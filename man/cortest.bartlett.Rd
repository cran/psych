\name{cortest.bartlett}
\alias{cortest.bartlett}

\title{Bartlett's test that a correlation matrix is an identity matrix }
\description{Bartlett (1951) proposed that -ln(det(R)*(N-1 - (2p+5)/6) was distributed as chi square if R were an identity matrix.  A useful test that residuals correlations are all zero.
}
\usage{
cortest.bartlett(R, n = NULL)
}

\arguments{
  \item{R}{A correlation matrix. (If R is not square, correlations are found and a warning is issued. }
  \item{n}{Sample size (if not specified, 100 is assumed.}
}
\details{More useful for pedagogical purposes than actual applications. The Bartlett test is asymptotically chi square distributed.
}
\value{
  \item{chisq}{Assymptotically chisquare}
  \item{p.value }{Of chi square}
  \item{df}{The degrees of freedom}

}
\references{ 
Bartlett, M. S., (1951), The Effect of Standardization on a chi square Approximation in Factor Analysis, Biometrika, 38, 337-344.
	
 }
\author{William Revelle}

\seealso{ \code{\link{cortest.mat}}, \code{\link{cortest.normal}}, \code{\link{cortest.jennrich}}}
\examples{
set.seed(42)   
x <- matrix(rnorm(1000),ncol=10)
r <- cor(x)
cortest.bartlett(r)      #random data don't differ from an identity matrix
data(bfi)
cortest.bartlett(bfi)    #not an identity matrix
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
