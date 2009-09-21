\name{dia.rect}
\Rdversion{1.1}
\alias{dia.rect}
\alias{dia.ellipse}
\alias{dia.ellipse1}
\alias{dia.arrow}
\alias{dia.curve}
\alias{dia.self}

\title{Helper functions for drawing path model diagrams}
\description{Path models are used to describe structural equation models or cluster analytic output.  These functions provide the primitives for drawing path models.  Used as a substitute for some of the functionality of Rgraphviz.}

\usage{
dia.rect(x, y = NULL, labels = NULL,  cex = 1,  xlim = c(0, 1), ylim = c(0, 1), ...)
 dia.ellipse(x, y = NULL, labels = NULL, cex=1,e.size=.05, xlim=c(0,1), ylim=c(0,1),  ...) 
dia.ellipse1(x,y,segments=51,e.size=.05,xlim=c(0,1),ylim=c(0,1),...)
dia.arrow(from,to,labels=NULL,scale=1,radius1=0,radius2=0,cex=1,...)
dia.curve(from,to,labels=NULL,scale=1,...)
dia.self(location,labels=NULL,scale=.8,side=2,...)
}

\arguments{
  \item{x}{x coordinate of a rectangle or ellipse}
  \item{y}{y coordinate of a rectangle or ellipse}
  \item{radius1}{start  the arrows at the endpoint - radius1 }
  \item{radius2}{stop the arrows at the endpoint - radius2 }
  \item{e.size}{The size of the ellipse (scaled by the number of variables}
  \item{labels}{Text to insert in rectangle, ellipse, or arrow}
  \item{cex}{adjust the text size}
  \item{segments}{precision with which to drawn the ellipses}
  \item{scale}{modifies size of rectangle and ellipse}
  \item{from}{arrows and curves go from }
  \item{to}{arrows and curves go to}
  \item{location}{where is the rectangle?}
  \item{xlim}{default ranges}
  \item{ylim}{default ranges}
  \item{side}{Which side of boxes should errors appear}
  \item{\dots}{Most graphic parameters may be passed here}

}
\details{These functions are the graphic primitives used by \code{\link{fa.diagram}}, \code{\link{omega.diagram}}, and \code{\link{ICLUST.diagram}}.  They create rectangles or ellipses surrounding text, connect them to straight or curved arrows, and can draw an arrow from and to the same rectangle.  

Each elipse and rectangle has a left, right, top and bottom coordinate that is used to connect the arrows. 

These functions were developed to get around the infelicities associated with trying to install Rgraphviz and graphviz.  

Better documentation will be added as these functions get improved.  Currently the functions are just a work around for Rgraphviz.
}
\value{Graphic output}

\author{William Revelle
}

\seealso{\code{\link{fa.diagram}}, \code{\link{omega.diagram}}, and \code{\link{ICLUST.diagram}}.
}
\examples{
v9 <- sim.hierarchical()
f3 <- fa(v9,3,rotate="cluster")
fa.diagram(f3,error=TRUE,side=3) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{hplot }% __ONLY ONE__ keyword per line
