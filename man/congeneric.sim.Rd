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
\details{When constructing examples for reliability analysis, it is convenient to simulate congeneric data structures.  These are the most simple of item structures, having just one factor.

The implied covariance matrix is just pattern \%*\% t(pattern). 
}
\value{
  \item{observed}{a matrix of test scores for n tests}
  \item{pattern }{The pattern matrix implied by the loadings and error variances}


}
\references{}
\author{ William Revelle }

}
\seealso{  \code{\link{item.sim}} }
\examples{

#test <- congeneric.sim()


## The function is currently defined as

function(N = 1000, loads = c(0.8, 0.7, 0.6, 0.5),err=NULL, short=TRUE) { 
 n <- length(loads) 
loading <- matrix(loads, nrow = n) 
error <- diag(1, nrow = n) 
if (!is.null(err)) {diag(error) <- err} else {
 diag(error) <- sqrt(1 - loading^2) }
pattern <- cbind(loading, error) 
colnames(pattern) <- c("theta", paste("e", seq(1:n), sep = "")) 
rownames(pattern) <- c(paste("V", seq(1:n), sep = "")) 
latent <- matrix(rnorm(N * (n + 1)), ncol = (n + 1)) 
observed <- latent %*% t(pattern) 
if (short) {return(observed)}  else result <- list(observed=observed,pattern=pattern)
 return(result)
} 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
