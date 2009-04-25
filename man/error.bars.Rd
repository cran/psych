\name{error.bars}
\alias{error.bars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot means and confidence intervals}
\description{One of the many functions in R to plot means and confidence intervals. Can be done using barplots if desired.  Can also be combined with such functions as boxplot to summarize distributions.  Means and standard errors are calculated from the raw data using \code{\link{describe}}.
}
\usage{
error.bars(x,stats=NULL, ylab = "Dependent Variable",xlab="Independent Variable", main=NULL, ylim = NULL, alpha=.05, labels = NULL, pos = NULL, arrow.len = 0.05, add = FALSE,bars=FALSE, ...)
}

\arguments{
  \item{x}{ A data frame or matrix of raw data }
  \item{stats}{Alternatively, a data.frame of descriptive stats from (e.g., describe)}
  \item{ylab}{y label}
  \item{xlab}{x label}
  \item{main}{title for figure}
  \item{ylim}{if specified, the limits for the plot, otherwise based upon the data}
  \item{alpha}{alpha level of confidence interval  -- defaults to 95\% confidence interval}
  \item{labels}{ X axis label }
  \item{pos}{where to place text: below, left, above, right}
  \item{arrow.len}{ How long should the top of the error bars be?}
  \item{add}{ add=FALSE, new plot, add=TRUE, just points and error bars}
  \item{bars}{bars=TRUE will draw a bar graph if you really want to do that}
  \item{\dots}{other parameters to pass to the plot function, e.g., typ="b" to draw lines, lty="dashed" to draw dashed lines}
}

\details{Drawing the mean +/- a confidence interval is a frequently used function when reporting experimental results. By default, the confidence interval is 1.96 standard errors.  
  
}
\value{Graphic output showing the means + x% confidence intervals.  For ci=1.96, and normal data, this will be the 95\% confidence region.  For ci=1, the 68\% confidence region.
 
}

\author{William Revelle}

\seealso{ \code{\link{error.crosses}} for two way error bars, \code{\link{error.bars.by}} for error bars for different groups }
\examples{
x <- replicate(20,rnorm(50))
boxplot(x,notch=TRUE,main="Notched boxplot with error bars")
error.bars(x,add=TRUE)
abline(h=0)

error.bars(attitude,alpha=.5,main="50 percent confidence limits") #another example
error.bars(attitude,bar=TRUE)  #show the use of bar graphs


#combine with a strip chart
stripchart(attitude,vertical=TRUE,method="jitter",main="Stripchart with 95 percent confidence limits")
error.bars(attitude,add=TRUE,arrow.len=.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
