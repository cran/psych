\name{draw.tetra}
\alias{draw.tetra}

\title{Draw a correlation ellipse and two normal curves to demonstrate tetrachoric correlation}

\description{A graphic of a correlation ellipse divided into 4 regions based upon x and y cutpoints on two normal distributions.  This is also an example of using the layout function.
}
\usage{
draw.tetra(r, t1, t2,shade=TRUE)
}

\arguments{
  \item{r}{the underlying Pearson correlation defines the shape of the ellipse}
  \item{t1}{X is cut at tau}
  \item{t2}{Y is cut at Tau}
  \item{shade}{shade the diagram (default is TRUE)}

}
\details{
A graphic demonstration of the \code{\link{tetrachoric}} correlation. Used for teaching purposes.  The default values are for a correlation of .5 with cuts at 1 and 1. Any other values are possible.  The code is also a demonstration of how to use the \code{\link{layout}} function for complex graphics using base graphics. 
}

\author{
William Revelle
}

\seealso{\code{\link{tetrachoric}}, \code{\link{irt.fa}}
}
\examples{
if(require(mvtnorm)) {
draw.tetra(.5,1,1)
draw.tetra(.8,2,1)} else {print("draw.tetra requires the mvtnorm package")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{hplot}