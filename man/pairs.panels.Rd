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
pairs.panels(x, y, smooth = TRUE, scale = FALSE, digits = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame or matrix}
  \item{y}{an optional data.frame or matrix }
  \item{smooth}{TRUE draws loess smooths }
  \item{scale}{ TRUE scales the correlation font by the size of the absolute correlation. }
  \item{digits}{ the number of digits to show}
  \item{\dots}{other options for pairs }
}
\details{Shamelessly adapted from the pairs help page.  Uses panel.cor, panel.cor.scale, and panel.hist, all taken from the help pages for pairs.  
}
\value{a scatter plot matrix (SPLOM) is drawn in the graphic window. The lower off diagonal draws scatter plots, the diagonal histograms, the upper off diagonal reports the Pearson correlation (with pairwise deletion).
}
\author{ 
 \email{  revelle@northwestern.edu } 
 }
\seealso{  \code{\link{pairs}} }
\examples{

#pairs.panels(attitude)   #see the graphics window

## The function is currently defined as
function (x, y, smooth = TRUE, scale = FALSE, digits = 2, ...) 
{
    if (smooth) {
        if (scale) {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale, 
                lower.panel = panel.smooth, ...)
        }
        else {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor, 
                lower.panel = panel.smooth, ...)
        }
    }
    else {
        if (scale) {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale, 
                ...)
        }
        else {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor, 
                ...)
        }
    }
  }
  
  #
  ## The function is currently defined as
"panel.cor" 
function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex )
     }
     
 ## The function is currently defined as
 "panel.cor.scale" 
function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex * abs(r))
     }
 "panel.hist"    
 function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
  }
}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
