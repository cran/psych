\name{paired.r}
\alias{paired.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test the difference between paired correlations }
\description{   Test the difference between paired correlations. Given 3 variables, x, y, z,  is the correlation between xy different than that between xz?  If y and z are independent, this is a simple t-test.  But, if they are dependent, it is a bit more complicated. 
}
\usage{
paired.r(xy, xz, yz, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{r(xy) }
  \item{xz}{r(xz) }
  \item{yz}{r(yz) }
  \item{n}{Number of subjects }
}
\details{Find a t2 statistic for the difference of two dependent correlations.
}
\value{t2
}
\references{ }
\author{ William Revelle \cr
         Northwestern University \cr
         Evanston, Illinois \cr
 \email{  revelle@northwestern.edu           } \cr
 \url{http://personality-project.org/revelle.html}}

\examples{
paired.r(.5,.3, .4, 100)
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
