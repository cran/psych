\name{pairs.panels}
\alias{pairs.panels}
\alias{panel.cor}
\alias{panel.cor.scale}
\alias{panel.hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SPLOM, histograms and correlations for a data matrix}
\description{Adapted from the  help page for pairs, pairs.panels shows a scatter plot of matrices (SPLOM), with bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal. Useful for descriptive statistics of small data sets.
}
\usage{
pairs.panels(x, smooth = TRUE, scale = FALSE, digits = 2, pch=20, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix}
  \item{smooth}{TRUE draws loess smooths }
  \item{scale}{ TRUE scales the correlation font by the size of the absolute correlation. }
  \item{digits}{ the number of digits to show}
  \item{pch}{The plot character (defaults to 20 which is a '.').}
  \item{\dots}{other options for pairs }
}
\details{Shamelessly adapted from the pairs help page.  Uses panel.cor, panel.cor.scale, and panel.hist, all taken from the help pages for pairs.  
}
\value{a scatter plot matrix (SPLOM) is drawn in the graphic window. The lower off diagonal draws scatter plots, the diagonal histograms, the upper off diagonal reports the Pearson correlation (with pairwise deletion).
}

\seealso{  \code{\link{pairs}} }
\examples{

pairs.panels(attitude)   #see the graphics window


}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
