\name{error.crosses}
\alias{error.crosses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot x and y error bars }
\description{Given two vectors of data, plot the means and show standard errors in both X and Y directions. 
}
\usage{
error.crosses(x, y, labels = NULL, pos = NULL, arrow.len = 0.2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of summary statistics (from Describe) }
  \item{y}{ A second vector of summary statistics (also from Describe)}
  \item{labels}{ name the pair }
  \item{pos}{Labels are located where with respect to the mean?  }
  \item{arrow.len}{ Arrow length }
  \item{\dots}{ Other parameters for plot  }
}
\details{For an example of two way error bars describing the effects of mood manipulations upon positive and negative affect, see \url{http://personality-project.org/revelle/publications/happy-sad-appendix/FIG.A-6.pdf}

The second example shows how error crosses can be done for multiple variables where the grouping variable is found dynamically. 
}

\author{ William Revelle \cr     
 \email{revelle@northwestern.edu}
}

\seealso{To draw error bars for single variables  \code{\link{error.bars}}, or by groups \code{\link{error.bars.by}}, or to find descriptive statistics \code{\link{describe}}  or descriptive statistics by a grouping variable \code{\link{describe.by}}  }
\examples{

desc <- describe(attitude)
x <- desc[1,]
y <- desc[2,]
plot(x$mean,y$mean,xlab=rownames(x),ylab=rownames(y))   #in graphics window
error.crosses(x,y)    #in graphics window
#now for a bit more complicated plotting 
desc <- describe.by(attitude,(attitude[,7]>41)) #select a high and low group
g1 <- desc$'FALSE'
g2 <- desc$'TRUE'
plot(g1$mean,g2$mean,xlab = "Low Advance",ylab="High Advance",xlim=c(30,80),ylim=c(50,80))
error.crosses(g1,g2,labels=rownames(g1),pos=rep(1,7))
title("Attitudes grouped by high and low scores on Advance")

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
