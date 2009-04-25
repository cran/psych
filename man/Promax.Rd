\name{Promax}
\alias{Promax}
\alias{target.rot}
\title{ Perform promax or targeted rotations and return the inter factor angles}
\description{promax is an oblique rotation function introduced by Hendrickson and White (1964) and implemented in the promax function in the stats package.  Unfortunately, promax does not report the inter factor correlations.  Promax does.  target.rot does general target rotations to an arbitrary target matrix. The default target rotation is for an independent cluster solution.
}
\usage{
Promax(x, m = 4)
target.rot(x,keys=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A loadings matrix}
  \item{m}{the power to which to raise the varimax loadings (for Promax)}
  \item{keys}{An arbitrary target matrix, can be composed of  any weights, but probably -1,0, 1 weights.  If missing, the target is the independent cluster structure determined by assigning every item to it's highest loaded factor.}
 
  
}
\details{Promax is a very direct adaptation of the stats::promax function.  The addition is that it will return the interfactor correlations as well as the loadings and rotation matrix. 

In addition, it will take output from either the factanal, \code{\link{factor.pa}}, \code{\link{factor.minres}} or \code{\link{principal}}  functions and select just the loadings matrix for analysis.

The target.rot function is an adaptation of a function of Michael Browne's to do rotations to arbitrary target matrices.  Suggested by Pat Shrout. 

The default for target.rot is to rotate to an independent cluster structure (every items is assigned to a group with its highest loading.) 

target.rot will not handle targets that have linear dependencies (e.g., a pure bifactor model where there is a g loading and a group factor for all variables).

}
\value{
  \item{loadings }{Oblique factor loadings}
  \item{rotmat}{The rotation matrix applied to the original loadings to produce the promax soluion or the targeted matrix}
  \item{Phi}{The interfactor correlation matrix}
}
\references{Hendrickson, A. E. and  White, P. O, 1964, British Journal of Statistical Psychology, 17, 65-70.
}
\author{William Revelle }
\note{A direct adaptation of the stats:promax function following suggestions to the R-help list by Ulrich Keller and John Fox. Furthe  modified to do  targeted rotation similar to a function of Michael Browne.

 
}
\seealso{ \code{\link{promax}}, \code{\link{factor.pa}}, \code{\link{factor.minres}}, or \code{\link{principal}} }
\examples{
jen <- sim.hierarchical()
f3 <- factor.minres(jen,3)
Promax(f3)
target.rot(f3)
m3 <- factanal(covmat=jen,factors=3)
Promax(m3)  #example of taking the output from factanal
#compare this rotation with the solution from a targeted rotation aimed for an independent cluster solution
target.rot(m3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models }% __ONLY ONE__ keyword per line
