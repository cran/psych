\name{error.bars}
\alias{error.bars}
\title{Plot means and confidence intervals}
\description{One of the many functions in R to plot means and confidence intervals. Can be done using barplots if desired.  Can also be combined with such functions as boxplot to summarize distributions.  Means and standard errors are calculated from the raw data using \code{\link{describe}}. Alternatively, plots of means +/- one standard deviation may be drawn.
}
\usage{
error.bars(x,stats=NULL, ylab = "Dependent Variable",xlab="Independent Variable",
  main=NULL, ylim = NULL, alpha=.05,sd=FALSE, labels = NULL, pos = NULL, 
    arrow.len = 0.05,arrow.col="black", add = FALSE,bars=FALSE,within=FALSE, ...)
}

\arguments{
  \item{x}{ A data frame or matrix of raw data }
  \item{stats}{Alternatively, a data.frame of descriptive stats from (e.g., describe)}
  \item{ylab}{y label}
  \item{xlab}{x label}
  \item{main}{title for figure}
  \item{ylim}{if specified, the limits for the plot, otherwise based upon the data}
  \item{alpha}{alpha level of confidence interval  -- defaults to 95\% confidence interval}
  \item{sd}{if TRUE, draw one  standard deviation instead of standard errors at the alpha level}
  \item{labels}{ X axis label }
  \item{pos}{where to place text: below, left, above, right}
  \item{arrow.len}{ How long should the top of the error bars be?}
  \item{arrow.col}{What color should the error bars be?}
  \item{add}{ add=FALSE, new plot, add=TRUE, just points and error bars}
  \item{bars}{bars=TRUE will draw a bar graph if you really want to do that}
  \item{within}{should the error variance of a variable be corrected by 1-SMC?}
  \item{\dots}{other parameters to pass to the plot function, e.g., typ="b" to draw lines, lty="dashed" to draw dashed lines}
}

\details{Drawing the mean +/- a confidence interval is a frequently used function when reporting experimental results. By default, the confidence interval is 1.96 standard errors of the t-distribution.  

If within=TRUE, the error bars are corrected for the correlation with the other variables by reducing the variance by a factor of (1-smc).  This allows for comparisons between variables.

The error bars are normally calculated from the data using the describe function.  If, alternatively, a matrix of statistics is provided with column headings of values, means, and se, then those values will be used for the plot (using the stats option).  However, in this case, the error bars will be one s.e. rather than a function of the alpha level.

If sd is TRUE, then the error bars will represent one standard deviation from the mean rather than be a function of alpha and the standard errors.
  
}
\value{Graphic output showing the means + x% confidence intervals.  For ci=1.96  and normal data, this will be the 95\% confidence region.  For ci=1, the 68\% confidence region. 

These confidence regions are based upon normal theory and do not take into account any skew in the variables.  More accurate confidence intervals could be found by resampling.

 
}

\author{William Revelle}

\seealso{ \code{\link{error.crosses}} for two way error bars, \code{\link{error.bars.by}} for error bars for different groups

In addition, as pointed out by Jim Lemon on the R-help news group, error bars or confidence intervals may be drawn using 
\tabular{ll}{
 function   \tab package \cr
 bar.err  \tab (agricolae) \cr
 plotCI \tab (gplots) \cr
 xYplot \tab(Hmisc) \cr
 dispersion \tab(plotrix) \cr
 plotCI \tab(plotrix) \cr
 }

For advice why not to draw bar graphs  with error bars, see  \url{http://biostat.mc.vanderbilt.edu/wiki/Main/DynamitePlots}

}
\examples{
x <- replicate(20,rnorm(50))
boxplot(x,notch=TRUE,main="Notched boxplot with error bars")
error.bars(x,add=TRUE)
abline(h=0)

error.bars(attitude,alpha=.5,main="50 percent confidence limits") #another example
error.bars(attitude,bar=TRUE)  #show the use of bar graphs


#combine with a strip chart and boxplot
stripchart(attitude,vertical=TRUE,method="jitter",jitter=.1,pch=19,
           main="Stripchart with 95 percent confidence limits")
boxplot(attitude,add=TRUE)
error.bars(attitude,add=TRUE,arrow.len=.2)

#use statistics from somewhere else
my.stats <- data.frame(values=c(1,4,8),mean=c(10,12,18),se=c(2,3,5))
error.bars(stats=my.stats,type="b",main="data with confidence intervals")
#note that in this case, the error bars are 1 s.e.  To modify that, change the s.e.

#Consider the case where we get stats from describe
temp <- describe(attitude)
error.bars(stats=temp)
#these error bars will be just one s.e.

#adjust the s.e. to vary by alpha level
alpha <- .05
temp[,"se"] <- temp[,"se"] * qt(1-alpha/2,temp[,"n"])
error.bars(stats=temp)
#show these do not differ from the other way by overlaying the two
error.bars(attitude,add=TRUE)


}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
