\name{error.bars.by}
\alias{error.bars.by}

\title{ Plot means and confidence intervals for multiple groups}
\description{One of the many functions in R to plot means and confidence intervals.  Meant mainly for demonstration purposes for showing the probabilty of replication from multiple samples.  Can also be combined with such functions as boxplot to summarize distributions.  Means and standard errors for each group are calculated using \code{\link{describe.by}}.
}
\usage{
error.bars.by(x,group,by.var=FALSE,x.cat=TRUE,ylab =NULL,xlab=NULL,main=NULL,ylim= NULL, alpha=.05,sd=FALSE, labels=NULL,pos=NULL,arrow.len=.05,add=FALSE,bars=FALSE,within=FALSE,colors=c("black","blue","red"),lty=NULL,...)
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
  \item{alpha}{alpha level of confidence interval.  Default is 1- alpha =95\% confidence interval}
  \item{sd}{sd=TRUE  will plot Standard Deviations instead of standard errors}
  \item{labels}{ X axis label }
  \item{pos}{where to place text: below, left, above, right}
  \item{arrow.len}{ How long should the top of the error bars be?}
  \item{add}{ add=FALSE, new plot, add=TRUE, just points and error bars}
  \item{bars}{Draw a barplot with error bars rather than a simple plot of the means}
  \item{within}{Should the s.e. be corrected by the correlation with the other variables?}
  \item{colors}{groups will be plotted in different colors (mod n.groups) }
  \item{lty}{line type may be specified in the case of not plotting by variables}
  \item{\dots}{other parameters to pass to the plot function, e.g., typ="b" to draw lines, lty="dashed" to draw dashed lines}

}
\details{Drawing the mean +/- a confidence interval is a frequently used function when reporting experimental results. By default, the confidence interval is 1.96 standard errors.  

This function was originally just a wrapper for \code{\link{error.bars}} but has been written to allow groups to be organized either as the x axis or as separate lines.

If desired, a barplot with error bars can be shown. Many find this type of plot to be uninformative (e.g.,   http://biostat.mc.vanderbilt.edu/DynamitePlots )
}
\value{Graphic output showing the means + x\% confidence intervals for each group.  For ci=1.96, and normal data, this will be the 95\% confidence region.  For ci=1, the 68\% confidence region.

These confidence regions are based upon normal theory and do not take into account any skew in the variables.  More accurate confidence intervals could be found by resampling.
}

\seealso{ See Also as \code{\link{error.crosses}}, \code{\link{error.bars}} }
\examples{

data(sat.act)
error.bars.by(sat.act[1:4],sat.act$gender)
error.bars.by(sat.act[5:6],sat.act$gender,bars=TRUE,labels=c("male","female"),main="SAT V and SAT Q by gender")  #draw a barplot

error.bars.by(sat.act[5:6],sat.act$education,bars=TRUE,xlab="Education",main="95 percent confidence limits of Sat V and Sat Q") 

error.bars.by(sat.act[5:6],sat.act$education,TRUE, xlab="Education")  #plot SAT V and SAT Q by education


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ hplot }% __ONLY ONE__ keyword per line
