\name{error.bars}
\alias{error.bars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot means and confidence intervals}
\description{One of the many functions in R to plot means and confidence intervals.  Meant mainly for demonstration purposes for showing the probabilty of replication from multiple samples.  Can also be combined with such functions as boxplot to summarize distributions.
}
\usage{
error.bars(x, ylab = "Dependent Variable",xlab="", ylim = NULL, ci = 1.96, labels = NULL, pos = NULL, arrow.len = 0.05, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame or matrix }
  \item{ylab}{y label}
  \item{xlab}{x label}
  \item{ylim}{if specified, the limits for the plot, otherwise based upon the data}
  \item{ci}{ What size confidence interval to use}
  \item{labels}{ X axis label }
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{arrow.len}{ How long should the top of the error bars be?}
  \item{add}{ add=FALSE, new plot, add=TRUE, just points and error bars}
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  
}
\value{
 
}

\author{William Revelle}

\seealso{  See Also as \code{\link{error.crosses}}, ~~~ }
\examples{
x <- matrix(rnorm(500),ncol=20)
error.bars(x)
#now do a boxplot and then add error bars
x.df <- as.data.frame(x)
boxplot(x.df)
error.bars(x.df, add=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
