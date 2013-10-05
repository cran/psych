\name{scatter.hist}
\alias{scatter.hist}

\title{Draw a scatter plot with associated X and Y histograms, densitie and correlation}
\description{Draw a X Y scatter plot with associated X and Y histograms with estimated densities.  Partly a demonstration of the use of layout. Also includes lowess smooth or linear model slope, as well as correlation.  Adapted from addicted to R example 78 
}
\usage{
scatter.hist(x,y=NULL,smooth=TRUE,ab=FALSE,correl=TRUE,density=TRUE,ellipse=TRUE,
  digits=2, method,cex.cor=1,title="Scatter plot + histograms",xlab=NULL,ylab=NULL,...)
}
\arguments{
  \item{x}{The X vector, or the first column of a  data.frame or matrix.  }
  \item{y}{The Y vector, of if X is a data.frame or matrix, the second column of X}
  \item{smooth}{if TRUE, then loess smooth it}
  \item{ab}{if TRUE, then show the best fitting linear fit}
  \item{correl}{TRUE: Show the correlation}
   \item{density}{TRUE: Show the estimated densities}
  \item{ellipse}{TRUE: draw 1 and 2 sigma ellipses and smooth}
  \item{digits}{How many digits to use if showing the correlation}
  \item{method}{Which method to use for correlation ("pearson","spearman","kendall") defaults to "pearson"}
  \item{cex.cor}{Adjustment for the size of the correlation}
  \item{xlab}{Label for the x axis}
  \item{ylab}{Label for the y axis} 
  \item{title}{An optional title}
  \item{\dots}{Other parameters for graphics}
}
\details{Just a straightforward application of layout and barplot, with some tricks taken from \code{\link{pairs.panels}}.  The various options allow for correlation ellipses (1 and 2 sigma from the mean), lowess smooths, linear fits, density curves on the histograms, and the value of the correlation.  ellipse = TRUE implies smooth = TRUE )
}
\author{William Revelle}
\note{
Adapted from Addicted to R example 78
}

\seealso{\code{\link{pairs.panels}} for multiple plots, \code{\link{multi.hist}} for multiple histograms.
}
\examples{
data(sat.act)
with(sat.act,scatter.hist(SATV,SATQ))
#or for something a bit more splashy
scatter.hist(sat.act[5:6],pch=(19+sat.act$gender),col=c("blue","red")[sat.act$gender])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ hplot }
