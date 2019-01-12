\name{ellipses}
\alias{ellipses}
\alias{minkowski}
\title{Plot data and 1 and 2 sigma correlation ellipses}
\description{For teaching correlation, it is useful to draw ellipses around the mean to reflect the correlation.  This variation of the ellipse function from John Fox's car package does so.  Input may be either two vectors or a matrix or data.frame.  In the latter cases, if the number of variables >2, then the ellipses are done in the \code{\link{pairs.panels}} function. Ellipses may be added to existing plots. 
The minkowski function is included as a generalized ellipse.
}
\usage{
ellipses(x, y = NULL, add = FALSE, smooth=TRUE, lm=FALSE,data=TRUE, n = 2,
 span=2/3, iter=3, col = "red", xlab =NULL,ylab= NULL,size=c(1,2), ...)
minkowski(r=2,add=FALSE,main=NULL,xl=1,yl=1) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector,matrix, or data.frame }
  \item{y}{Optional second vector }
  \item{add}{Should a new plot be created, or should it be added to?}
  \item{smooth}{smooth = TRUE -> draw a loess fit}
  \item{lm}{lm=TRUE -> draw the linear fit}
  \item{data}{data=TRUE implies draw the data points}
  \item{n}{Should 1 or 2 ellipses be drawn }
  \item{span}{averaging window parameter for the lowess fit}
  \item{iter}{iteration parameter for lowess}
  \item{col}{color of ellipses (default is red}
  \item{xlab}{label for the x axis}
  \item{ylab}{label for the y axis}
  \item{size}{The size of ellipses in sd units (defaults to 1 and 2)}
  \item{\dots}{Other parameters for plotting}
  \item{r}{r=1 draws a city block, r=2 is a Euclidean circle, r > 2 tends towards a square}
  \item{main}{title to use when drawing Minkowski circles}
  \item{xl}{stretch the x axis}
  \item{yl}{stretch the y axis}
  
}
\details{Ellipse dimensions are calculated from the correlation between the x and y variables and are scaled as sqrt(1+r) and sqrt(1-r). They are then scaled as size[1] and size[2] standard deviation units.   To scale  for 95 and 99 percent confidence use c(1.64,2.32)
}
\value{ A single plot (for 2 vectors or data frames with fewer than 3 variables.  Otherwise a call is made to \code{\link{pairs.panels}}. 
}
\references{  Galton, Francis (1888), Co-relations and their measurement. Proceedings of the Royal Society. London Series, 45, 135-145.
 }
\author{William Revelle }
\note{Adapted from John Fox's ellipse and data.ellipse functions.
}
\seealso{  \code{\link{pairs.panels}}  }
\examples{
data(galton)
ellipses(galton,lm=TRUE)
ellipses(galton$parent,galton$child,xlab="Mid Parent Height",
                  ylab="Child Height") #input are two vectors
data(sat.act)
ellipses(sat.act)  #shows the pairs.panels ellipses
minkowski(2,main="Minkowski circles")
minkowski(1,TRUE)
minkowski(4,TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line

