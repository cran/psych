\name{Promax}
\alias{Promax}
\title{ Perform promax rotation and return the inter factor angles}
\description{promax is an oblique rotation function introduced by Hendrickson and White (1964) and implemented in the promax function in the stats package.  Unfortunately, promax does not report the inter factor correlations.  Promax does. 
}
\usage{
Promax(x, m = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A loadings matrix}
  \item{m}{the power to which to raise the varimax loadings }
}
\details{This is a very direct adaptation of the stats::promax function.  The addition is that it will return the interfactor correlations as well as the loadings and rotation matrix. 

In addition, it will take output from either the factanal, factor.pa, or principal functions and select just the loadings matrix for analysis.
}
\value{
  \item{loadings }{Oblique factor loadings}
  \item{rotmat}{The rotation matrix applied to the original loadings to produce the promax soluion}
  \item{Phi}{The interfactor correlation matrix}
}
\references{Hendrickson, A. E. and  White, P. O, 1964, British Journal of Statistical Psychology, 17, 65-70.
}
\author{William Revelle }
\note{A direct adaptation of the stats:promax function following suggestions to the R-help list by Ulrich Keller and John Fox.

 
}
\seealso{ \code{\link{promax}}, \code{\link{factor.pa}} }
\examples{
jen <- sim.hierarchical()
f3 <- factor.pa(jen,3)
Promax(f3)
m3 <- factanal(covmat=jen,factors=3)
Promax(m3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models }% __ONLY ONE__ keyword per line
