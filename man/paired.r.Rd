\name{paired.r}
\alias{paired.r}
\title{ Test the difference between (un)paired correlations }
\description{   Test the difference between two (paired or unpaired) correlations. Given 3 variables, x, y, z,  is the correlation between xy different than that between xz?  If y and z are independent, this is a simple t-test of the z transformed rs.  But, if they are dependent, it is a bit more complicated. 
}
\usage{
paired.r(xy, xz, yz=NULL, n, n2=NULL,twotailed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{r(xy) }
  \item{xz}{r(xz) }
  \item{yz}{r(yz) }
  \item{n}{Number of subjects for first group}
  \item{n2}{Number of subjects in second group (if not equal to n)}
  \item{twotailed}{Calculate two or one tailed probability values}
  }
\details{To find the z of the difference between two independent correlations, first convert them to z scores using the Fisher r-z transform and then find the z of the difference between the two correlations.  The default assumption is that the group sizes are the same, but the test can be done for different size groups by specifying n2.

If the correlations are not independent (i.e., they are from the same sample) then the correlation with the third variable r(yz) must be specified. Find a t statistic for the difference of thee two dependent correlations.
}
\value{a list containing the calculated t or z values and the associated two (or one) tailed probability.
  \item{t}{t test of the difference between two dependent correlations}
  \item{p}{probability of the t or of the z}
  \item{z}{z test of the difference between two independent correlations}
}

\seealso{\code{\link{p.rep.r}}, \code{\link{cor.test}} } 
\author{ William Revelle \cr
         Northwestern University \cr
         Evanston, Illinois \cr
 \email{  revelle@northwestern.edu           } \cr
 \url{http://personality-project.org/revelle.html}}

\examples{
paired.r(.5,.3, .4, 100) #dependent correlations
paired.r(.5,.3,NULL,100) #independent correlations same sample size
paired.r(.5,.3,NULL, 100, 64) # independent correlations, different sample sizes
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
