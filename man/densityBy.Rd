\name{densityBy}
\alias{densityBy}
\alias{violinBy}

\title{Create a 'violin plot' or density plot of the distribution of a set of variables}

\description{Among the many ways to describe a data set, one is density plot for each value of a grouping variable and another is violin plot of multiple variables.  A density plot shows the density for different groups to show effect sizes. A violin plot is similar 
to a box plot but shows the actual distribution.
Median and 25th and 75th percentile lines are added to the display. If a grouping variable is specified, violinBy will draw violin plots for each variable and for each group. }
\usage{
violinBy(x,var=NULL,grp=NULL,grp.name=NULL,ylab="Observed",xlab="",main="Density plot",
alpha= 1,adjust=1,restrict=TRUE,xlim=NULL,add=FALSE,col=NULL,pch=20,scale=NULL,...) 
 
densityBy(x,var=NULL,grp=NULL,freq=FALSE,col=c("blue","red","black"),alpha=.5,adjust=1,
 xlab="Variable", ylab="Density",main="Density Plot")

}


\arguments{
  \item{x}{A matrix or data.frame}
  \item{var}{The variable(s) to display}
  \item{grp}{A grouping variable}
  \item{grp.name}{If the grouping variable is specified, then what names should be give to the group? Defaults to 1:ngrp}
  \item{ylab}{The y label}
  \item{xlab}{The x label}
  \item{main}{Figure title}
  \item{alpha}{A degree of transparency (0=transparent ... 1 not transparent)}
  \item{adjust}{Allows smoothing of density histograms when plotting variables like height}
  \item{freq}{if TRUE, then plot frequencies (n * density)}
  \item{restrict}{Restrict the density to the observed max and min of the data}
  \item{xlim}{if not specified, will be .5 beyond the number of variables}
  \item{add}{Allows overplotting}
  \item{col}{Allows for specification of colours.  The default for 2 groups is 
          blue and red, for more group levels, rainbows.}
 \item{pch}{The plot character for the mean is by default a small filled circle.  To not show the mean, use pch=NA }
 \item{scale}{If NULL, scale the widths by the square root of sample size, otherwise scale by the value supplied.}
  \item{\dots}{Other graphic parameters}
}
\details{
Describe the data using a violin plot. Change alpha to modify the shading.  
The grp variable may be used to draw separate violin plots for each of multiple groups.
}
\value{
The density (y axis) by value (x axis)  of the data  (for densityBy) or a violin plot for each variable (perhaps broken down by groups)
}
\author{William Revelle
}
\note{Nothing yet}
\seealso{  \code{\link{describe}}, \code{\link{describeBy}} and  \code{\link{statsBy}} for descriptive statistics and \code{\link{error.bars}} \code{\link{error.bars.by}} and  \code{\link{bi.bars}} for graphic displays }
\examples{
violinBy(psychTools::bfi[1:5])
violinBy(psychTools::bfi,var=1:5,grp ="gender",grp.name=c("M","F"))
densityBy(sat.act,"SATV","gender")
  
}
\keyword{ multivariate }
\keyword{ hplot }
