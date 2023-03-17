\name{sim.congeneric}
\alias{congeneric.sim}
\alias{sim.congeneric}
\alias{make.congeneric}
\title{ Simulate a congeneric data set with or without minor factors }
\description{Classical Test Theory (CTT) considers four or more tests to be congenerically equivalent if all tests may be expressed in terms of one factor and a residual error.  Parallel tests are the special case where (usually two) tests have equal factor loadings.  Tau equivalent tests have equal factor loadings but may have unequal errors.  Congeneric tests may differ in both factor loading and error variances.  Minor factors may be added as systematic but trivial disturbances 
}
\usage{
sim.congeneric(loads = c(0.8, 0.7, 0.6, 0.5),N = NULL,  err=NULL, short = TRUE, 
              categorical=FALSE, low=-3,high=3,cuts=NULL,minor=FALSE,fsmall = c(-.2,.2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{How many subjects to simulate. If NULL, return the population model }
  \item{loads}{ A vector of factor loadings for the tests  }
  \item{err}{A vector of error variances -- if NULL then error = 1 - loading 2}
  \item{short}{short=TRUE: Just give the test correlations, short=FALSE, report observed test scores as well as the implied pattern matrix}
  \item{categorical}{ continuous or categorical (discrete) variables.  }
  \item{low}{ values less than low are forced to low }
  \item{high}{ values greater than high are forced to high  }
  \item{cuts}{If specified, and categorical = TRUE, will cut the resulting continuous output at the value of cuts}
  \item{minor}{Should n/2 minor factors be added (see Maccallum and Tucker, 1991)}
  \item{fsmall}{nvar/2 small factors are generated with loadings sampled from fsmall e.g. (-.2,0,.2)}
}
\details{When constructing examples for reliability analysis, it is convenient to simulate congeneric data structures.  These are the most simple of item structures, having just one factor. Mainly used for a discussion of reliability theory as well as factor score estimates. 

Maccallum and Tucker (1991) suggest that factor models should include minor factors, that at not random error but unspecifed by the basic model.  This option has been added in November, 2022.

The implied covariance matrix is just pattern \%*\% t(pattern). 
}
\value{
  \item{model}{The implied population correlation matrix if N=NULL or short=FALSE, otherwise the sample correlation matrix}
  \item{pattern }{The pattern matrix implied by the loadings and error variances}
  \item{r}{The sample correlation matrix for long output}
  \item{observed}{a matrix of test scores for n tests}
  \item{latent}{The latent trait and error scores }


}
\references{Revelle, W. (in prep) An introduction to psychometric theory with applications in R. To be published by Springer.  (working draft available at  \url{https://personality-project.org/r/book/} 

MacCallum, R. C., & Tucker, L. R. (1991). Representing sources of error in the 
common-factormodel: Implications for theory and practice. Psychological Bulletin, 109(3), 502-511.}
\author{ William Revelle }

\seealso{  \code{\link{item.sim}} for other simulations, \code{\link{fa}} for an example of factor scores, \code{\link{irt.fa}} and \code{\link{polychoric}} for the treatment of item data with discrete values.}

\examples{
test <- sim.congeneric(c(.9,.8,.7,.6))   #just the population matrix
test <- sim.congeneric(c(.9,.8,.7,.6),N=100)   # a sample correlation matrix
test <- sim.congeneric(short=FALSE, N=100)
round(cor(test$observed),2) # show  a congeneric correlation matrix
f1=fa(test$observed,scores=TRUE)
round(cor(f1$scores,test$latent),2)  
     #factor score estimates are correlated with but not equal to the factor scores
set.seed(42)
#500 responses to 4 discrete items
items <- sim.congeneric(N=500,short=FALSE,low=-2,high=2,categorical=TRUE) 
d4 <- irt.fa(items$observed)  #item response analysis of congeneric measures



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
\keyword{datagen}
