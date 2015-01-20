\name{cortest.mat}
\alias{cortest.normal}
\alias{cortest.mat}
\alias{cortest.jennrich}
\alias{cortest}

\title{Chi square tests of whether a single matrix is an identity matrix, or a pair of matrices are equal. }
\description{Steiger (1980) pointed out that the sum of the squared elements of a correlation matrix, or the Fisher z score equivalents, is distributed as chi square under the null hypothesis that the values are zero (i.e., elements of the identity matrix).  This is particularly useful for examining whether correlations in a single matrix differ from zero or for comparing two matrices. Jennrich (1970) also examined tests of differences between matrices.
}
\usage{
cortest.normal(R1, R2 = NULL, n1 = NULL, n2 = NULL, fisher = TRUE) #the steiger test
cortest(R1,R2=NULL,n1=NULL,n2 = NULL, fisher = TRUE,cor=TRUE) #same as cortest.normal 
cortest.jennrich(R1,R2,n1=NULL, n2=NULL)  #the Jennrich test
cortest.mat(R1,R2=NULL,n1=NULL,n2 = NULL) #an alternative test

}

\arguments{
  \item{R1}{A correlation matrix. (If R1 is not rectangular, and cor=TRUE, the correlations are found). }
  \item{R2}{A correlation matrix.  If R2 is not rectangular, and cor=TRUE, the correlations are found. If R2 is NULL, then the test is just whether R1 is an identity matrix. }
  \item{n1}{Sample size of R1 }
  \item{n2}{Sample size of R2 }
  \item{fisher}{Fisher z transform the correlations? }
  \item{cor}{By default, if the input matrices are not symmetric, they are converted to correlation matrices.  That is, they are treated as if they were the raw data.  If cor=FALSE, then the input matrices are taken to be correlation matrices.} 
}
\details{There are several ways to test if a matrix is the identity matrix. The most well known is the chi square test  of Bartlett (1951) and Box (1949). A very straightforward test, discussed by Steiger (1980) is to find the sum of the squared correlations or the sum of the squared Fisher transformed correlations.  Under the null hypothesis that all the correlations are equal, this sum is distributed as chi square.  This is implemented in 
\code{\link{cortest}} and \code{\link{cortest.normal}}

Yet another test, is the Jennrich(1970) test of the equality of two matrices. This compares the differences between two matrices to the averages of two matrices using a chi square test. This is implemented in \code{\link{cortest.jennrich}}.

Yet another option \code{\link{cortest.mat}} is to compare the two matrices using an approach analogous to that used in evaluating the adequacy of a factor model.  In factor analysis, the maximum likelihood fit statistic is 
  \cr
  \eqn{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^{-1} R|) - n.items}{f = log(trace ((FF'+U2)^{-1} R) -  log(|(FF'+U2)^-1 R|) - n.items}. 
  
This in turn is converted to a chi square 

 \eqn{\chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f }{chi^2 = (n.obs - 1 - (2 * p + 5)/6 - (2 * factors)/3)) * f } (see \code{\link{fa}}.)
 
 That is, the model (M = FF' + U2) is compared to the original correlation matrix (R) by a function of \eqn{M^{-1} R}.  By analogy, in the case of two matrices, A and B, \code{\link{cortest.mat}} finds the chi squares associated with \eqn{A^{-1}B} and \eqn{A B^{-1}}.  The sum of these two \eqn{\chi^2} will also be a \eqn{\chi^2} but with twice the degrees of freedom.  
 
 
}
\value{
  \item{chi2}{The chi square statistic}
  \item{df}{Degrees of freedom for the Chi Square}
  \item{prob}{The probability of observing the Chi Square under the null hypothesis.}
  
}
\references{ Steiger, James H. (1980) Testing pattern hypotheses on correlation matrices: alternative statistics and some empirical results. Multivariate Behavioral Research, 15, 335-352.

Jennrich, Robert I. (1970) An Asymptotic \eqn{\chi^2} Test for the Equality of Two Correlation Matrices. Journal of the American Statistical Association, 65, 904-912.
}


\author{ William Revelle }
\note{ Both the  cortest.jennrich  and cortest.normal  are probably overly stringent.  The ChiSquare values for pairs of random samples from the same population are larger than would be expected.  This is a good test for rejecting the null of no differences. 
}

\seealso{\code{\link{cortest.bartlett}} }
\examples{
x <- matrix(rnorm(1000),ncol=10)
cortest.normal(x)  #just test if this matrix is an identity
x <- sim.congeneric(loads =c(.9,.8,.7,.6,.5),N=1000,short=FALSE)
y <- sim.congeneric(loads =c(.9,.8,.7,.6,.5),N=1000,short=FALSE)
cortest.normal(x$r,y$r,n1=1000,n2=1000) #The Steiger test
cortest.jennrich(x$r,y$r,n1=100,n2=1000) # The Jennrich test
cortest.mat(x$r,y$r,n1=1000,n2=1000)   #twice the degrees of freedom as the Jennrich


}

\keyword{ multivariate }
