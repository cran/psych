\name{pairs.panels}
\alias{pairs.panels}
\alias{panel.cor}
\alias{panel.cor.scale}
\alias{panel.hist}
\alias{panel.lm}
\alias{panel.lm.ellipse}
\alias{panel.hist.density}
\alias{panel.ellipse}
\alias{panel.smoother}

\title{SPLOM, histograms and correlations for a data matrix}
\description{Adapted from the  help page for pairs, pairs.panels shows a scatter plot of matrices (SPLOM), with bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal. Useful for descriptive statistics of small data sets.  If lm=TRUE, linear regression fits are shown for both y by x and x by y.  Correlation ellipses are also shown. Points may be given different colors depending upon some grouping variable.
}
\usage{
pairs.panels(x, smooth = TRUE, scale = FALSE, density=TRUE,ellipses=TRUE,digits = 2, pch = 20,lm=FALSE, jiggle=FALSE,factor=2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix}
  \item{smooth}{TRUE draws loess smooths }
  \item{scale}{ TRUE scales the correlation font by the size of the absolute correlation. }
  \item{density}{TRUE shows the density plots as well as histograms}
  \item{ellipses}{TRUE draws correlation ellipses}
  \item{lm}{Plot the linear fit rather than the LOESS smoothed fits.}
  \item{digits}{ the number of digits to show}
  \item{pch}{The plot character (defaults to 20 which is a '.').}
  \item{jiggle}{Should the points be jittered before plotting?}
  \item{factor}{factor for jittering (1-5)}
  \item{\dots}{other options for pairs }
}
\details{Shamelessly adapted from the pairs help page.  Uses panel.cor, panel.cor.scale, and panel.hist, all taken from the help pages for pairs. Also adapts the ellipse function from John Fox's car package. 

\code{\link{pairs.panels}} is most useful when the number of variables to plot is less than about 6-10. It is particularly useful for an initial overview of the data.

To show different groups with different colors, use a plot character (pch) between 21 and 25 and then set the background color to vary by group. (See the second example).

}
\value{A scatter plot matrix (SPLOM) is drawn in the graphic window. The lower off diagonal draws scatter plots, the diagonal histograms, the upper off diagonal reports the Pearson correlation (with pairwise deletion).

If lm=TRUE, then the scatter plots are drawn above and below the diagonal, each with a linear regression fit.  Useful to show the difference between regression lines.
}

\seealso{  \code{\link{pairs}} }
\examples{

pairs.panels(attitude)   #see the graphics window
data(iris)
pairs.panels(iris[1:4],bg=c("red","yellow","blue")[iris$Species],pch=21,main="Fisher Iris data by Species") #to show color grouping


}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
