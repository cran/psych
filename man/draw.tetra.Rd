\name{draw.tetra}
\alias{draw.tetra}
\alias{draw.cor}

\title{Draw a correlation ellipse and two normal curves to demonstrate tetrachoric correlation}

\description{A graphic of a correlation ellipse divided into 4 regions based upon x and y cutpoints on two normal distributions.  This is also an example of using the layout function. Draw a bivariate density plot to show how tetrachorics work.
}
\usage{
draw.tetra(r, t1, t2,shade=TRUE)
draw.cor(r=.5,expand=10,theta=30,phi=30,N=101,nbcol=30,box=TRUE,
main="Bivariate density  rho = ",cuts=NULL,all=TRUE,ellipses=TRUE,ze=.15)
}

\arguments{
  \item{r}{the underlying Pearson correlation defines the shape of the ellipse}
  \item{t1}{X is cut at tau}
  \item{t2}{Y is cut at Tau}
  \item{shade}{shade the diagram (default is TRUE)}
  \item{expand}{The relative height of the z axis}
  \item{theta}{The angle to rotate the x-y plane}
  \item{phi}{The angle above the plane to view the graph}
  \item{N}{The grid resolution}
  \item{nbcol}{The color resolution}
  \item{box}{Draw the axes}
  \item{main}{The main title}
  \item{cuts}{Should the graphic show cuts (e.g., cuts=c(0,0))}
  \item{all}{Show all four parts of the tetrachoric}
  \item{ellipses}{Draw a correlation ellipse}
  \item{ze}{height of the ellipse if requested}

}
\details{
A graphic demonstration of the \code{\link{tetrachoric}} correlation. Used for teaching purposes.  The default values are for a correlation of .5 with cuts at 1 and 1. Any other values are possible.  The code is also a demonstration of how to use the \code{\link{layout}} function for complex graphics using base graphics. 
}

\author{
William Revelle
}

\seealso{\code{\link{tetrachoric}} to find tetrachoric correlations, \code{\link{irt.fa}} and \code{\link{fa.poly}} to use them in factor analyses, \code{\link{scatter.hist}} to show correlations and histograms. }

\examples{
if(require(mvtnorm)) {
draw.tetra(.5,1,1)
draw.tetra(.8,2,1)} else {print("draw.tetra requires the mvtnorm package")
draw.cor(.5,cuts=c(0,0))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{hplot}