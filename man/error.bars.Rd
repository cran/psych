\name{error.bars}
\alias{error.bars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot means and confidence intervals}
\description{One of the many functions in R to plot means and confidence intervals.  Meant mainly for demonstration purposes for showing the probabilty of replication from multiple samples.  Can also be combined with such functions as boxplot to summarize distributions.  Means and standard errors are calculated from the raw data using \code{\link{describe}}.
}
\usage{
error.bars(x, ylab = "Dependent Variable",xlab="Independent Variable", main="95\% confidence limits", ylim = NULL, ci = 1.96, labels = NULL, pos = NULL, arrow.len = 0.05, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame or matrix }
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
  
}
\value{Graphic output showing the means + x% confidence intervals.  For ci=1.96, and normal data, this will be the 95\% confidence region.  For ci=1, the 68\% confidence region.
 
}

\author{William Revelle}

\seealso{  See Also as \code{\link{error.crosses}}, \code{\link{error.bars.by}} }
\examples{
x <- matrix(rnorm(500),ncol=20)
error.bars(x)
#now do a boxplot and then add error bars
x.df <- as.data.frame(x)
boxplot(x.df)
error.bars(x.df, add=TRUE)

error.bars(attitude) #another example

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line