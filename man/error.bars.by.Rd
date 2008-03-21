\name{error.bars.by}
\alias{error.bars.by}

\title{ Plot means and confidence intervals for multiple groups}
\description{One of the many functions in R to plot means and confidence intervals.  Meant mainly for demonstration purposes for showing the probabilty of replication from multiple samples.  Can also be combined with such functions as boxplot to summarize distributions.  Means and standard errors for each group are calculated using \code{\link{describe.by}}.
}
\usage{
error.bars.by(x, group,by.var = FALSE,x.cat=TRUE, ylab = "NULL", xlab = "NULL",main="95\% Confidence Intervals", ylim = NULL, ci = 1.96, labels = NULL, pos = NULL, arrow.len = 0.05, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame or matrix }
  \item{group}{A grouping variable}
  \item{by.var}{A different line for each group (default) or each variable}
  \item{x.cat}{Is the grouping variable categorical (TRUE) or continuous (FALSE}
  \item{ylab}{y label}
  \item{xlab}{x label}
  \item{main}{title for figure}
  \item{ylim}{if specified, the limits for the plot, otherwise based upon the data}
  \item{ci}{ What size confidence interval to use}
  \item{labels}{ X axis label }
  \item{pos}{where to place text: below, left, above, right}
  \item{arrow.len}{ How long should the top of the error bars be?}
  \item{add}{ add=FALSE, new plot, add=TRUE, just points and error bars}
  \item{\dots}{other parameters to pass to the plot function, e.g., typ="b" to draw lines, lty="dashed" to draw dashed lines}

}
\details{Drawing the mean +/- a confidence interval is a frequently used function when reporting experimental results. By default, the confidence interval is 1.96 standard errors.  

This function is a wrapper for \code{\link{error.bars}} and allows groups to be organized either as the x axis or as separate lines.
}
\value{Graphic output showing the means + x\% confidence intervals for each group.  For ci=1.96, and normal data, this will be the 95\% confidence region.  For ci=1, the 68\% confidence region.
}

\seealso{ See Also as \code{\link{error.crosses}}, \code{\link{error.bars}} }
\examples{

x <- matrix(rnorm(500),ncol=20)
y <- sample(4,25 ,replace=TRUE)
x <- x+y
error.bars.by(x,y)
error.bars.by(x,y,TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ hplot }% __ONLY ONE__ keyword per line
