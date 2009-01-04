\name{cortest.mat}
\alias{cortest.normal}
\alias{cortest.mat}
\alias{cortest.jennrich}
\alias{cortest}

\title{Chi square tests of whether a single matrix is an identity matrix, or a pair of matrices are equal. }
\description{Steiger (1980) pointed out that the sum of the squared elements of a correlation matrix, or the Fisher z score equivalents, is distributed as chi square under the null hypothesis that the values are zero (i.e., elements of the identity matrix).  This is particularly useful for examining whether correlations in a single matrix differ from zero or for comparing two matrices. Jennrich (1970) also examined tests of differences between matrices.
}
\usage{
cortest.normal(R1, R2 = NULL, n1 = NULL, n2 = NULL, fisher = TRUE)
cortest(R1,R2=NULL,n1=NULL,n2 = NULL, fisher = TRUE)   #same as cortest.normal 
cortest.mat(R1,R2=NULL,n1=NULL,n2 = NULL)
cortest.jennrich(R1,R2,n1=NULL, n2=NULL)
}

\arguments{
  \item{R1}{A correlation matrix. (If R1 is not rectangular, the correlations are found). }
  \item{R2}{A correlation matrix.  If R2 is not rectangular, the correlations are found. If R2 is NULL, then the test is just whether R1 is an identity matrix. }
  \item{n1}{Sample size of R1 }
  \item{n2}{Sample size of R2 }
  \item{fisher}{Fisher z transform the correlations? }
}
\details{There are several ways to test if a matrix is the identity matrix. The most well known is the chi square test  of Bartlett (1951) and Box (1949). A very straightforward test, discussed by Steiger (1980) is to find the sum of the squared correlations or the sum of the squared Fisher transformed correlations.  Under the null hypothesis that all the correlations are equal, this sum is distributed as chi square.   

Yet another test, is the Jennrich(1970) test of the equality of two matrices.
}
\value{
  \item{chi2}{The chi square statistic}
  \item{df}{Degrees of freedom for the Chi Square}
  \item{prob}{The probability of observing the Chi Square under the null hypothesis.}
  
}
\references{ Steiger, James H. (1980) Testing pattern hypotheses on correlation matrices: alternative statistics and some empirical results. Multivariate Behavioral Research, 15, 335-352.}
\author{ William Revelle }
\note{ Both the  cortest.jennrich  and cortest.normal  are probably overly stringent.  The ChiSquare values for pairs of random samples from the same population are larger than would be expected.  This is a good test for rejecting the null of no differences. 
}

\seealso{\code{\link{cortest.bartlett}} }
\examples{
x <- matrix(rnorm(1000),ncol=10)
y <- matrix(rnorm(500),ncol=10)
cortest.normal(x)  #just test if this matrix is an identity
cortest.normal(x,y) #do these two matrices differ?
cortest.mat(x)
cortest.mat(x,y)   #twice the degrees of freedom as the Jennrich
cortest.jennrich(x,y) #

}

\keyword{ multivariate }