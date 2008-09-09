\name{congeneric.sim}
\alias{congeneric.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a congeneric data set }
\description{Classical Test Theory (CTT) considers four or more tests to be congenerically equivalent if all tests may be expressed in terms of one factor and a residual error.  Parallel tests are the special case where (usually two) tests have equal factor loadings.  Tau equivalent tests have equal factor loadings but may have unequal errors.  Congeneric tests may differ in both factor loading and error variances.
}
\usage{
congeneric.sim(N = 1000, loads = c(0.8, 0.7, 0.6, 0.5), err=NULL, short = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{How many subjects to simulate }
  \item{loads}{ A vector of factor loadings for the tests  }
  \item{err}{A vector of error variances -- if NULL then error = 1 - loading 2}
  \item{short}{short=TRUE: Just give the test scores, short=FALSE, report observed test scores as well as the implied pattern matrix}
}
\details{When constructing examples for reliability analysis, it is convenient to simulate congeneric data structures.  These are the most simple of item structures, having just one factor. Mainly used for a discussion of reliability theory as well as factor score estimates.

The implied covariance matrix is just pattern \%*\% t(pattern). 
}
\value{
  \item{model}{The implied population correlation matrix}
  \item{pattern }{The pattern matrix implied by the loadings and error variances}
  \item{observed}{a matrix of test scores for n tests}
  \item{latent}{The latent trait and error scores }


}
\references{}
\author{ William Revelle }

\seealso{  \code{\link{item.sim}} for other simulations, \code{\link{factor.pa}} for an example of factor scores.}
\examples{

test <- congeneric.sim(short=FALSE)
round(cor(test$observed),2) # show  a congeneric correlation matrix
f1=factor.pa(test$observed,1,scores=TRUE)
round(cor(f1$scores,test$latent),2)  #factor score estimates are correlated with but not equal to the factor scores



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
\keyword{datagen}
