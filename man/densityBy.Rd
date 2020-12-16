\name{densityBy}
\alias{densityBy}
\alias{violinBy}
\alias{violin}

\title{Create a 'violin plot' or density plot of the distribution of a set of variables}

\description{Among the many ways to describe a data set, one is a density plot for each value of a grouping variable and another is violin plot of multiple variables.  A density plot shows the density for different groups to show effect sizes. A violin plot is similar 
to a box plot but shows the actual distribution.
Median and 25th and 75th percentile lines are added to the display. If a grouping variable is specified, violinBy will draw violin plots for each variable and for each group. Data points may be drawn as well.
 }
\usage{

violin(x,data=NULL,var=NULL,grp=NULL,grp.name=NULL,xlab=NULL,ylab=NULL,
main="Density plot",vertical=TRUE,dots=FALSE,jitter=.05,alpha=1,
errors=FALSE,eyes=TRUE,adjust=1, restrict=TRUE,xlim=NULL, add=FALSE, 
col=NULL,pch=20,scale=NULL,...) 
 
violinBy(x,var=NULL,grp=NULL,data=NULL,grp.name=NULL,xlab=NULL,ylab=NULL,
main="Density plot",vertical=TRUE,dots=FALSE,jitter=.05,alpha= 1,
errors=FALSE,eyes=TRUE,adjust=1, restrict=TRUE,xlim=NULL, add=FALSE, 
col=NULL,pch=20,scale=NULL,...) 
 
densityBy(x,var=NULL,grp=NULL,data=NULL,freq=FALSE,col=c("blue","red","black"),alpha=.5, 
 adjust=1,ylim=NULL,xlim=NULL,xlab="Variable", ylab="Density",
       main="Density Plot",legend=NULL)

}


\arguments{
  \item{x}{A matrix or data.frame (can be expressed in formula input)}
  \item{var}{The variable(s) to display}
  \item{grp}{The grouping variable(s) }
  \item{data}{The name of the data object if using formula input}
  \item{grp.name}{If the grouping variable is specified, then what names should be give to the group? Defaults to 1:ngrp}
  \item{ylab}{The y label}
  \item{xlab}{The x label}
  \item{main}{Figure title}
  \item{vertical}{If TRUE, plot the violins vertically, otherwise, horizontonally}
  \item{dots}{if TRUE, add a stripchart with the data points}
  \item{jitter}{If doing a stripchart, then jitter the points this much}
  \item{errors}{If TRUE, add error bars or cats eyes to the violins}\
  \item{eyes}{if TRUE and errors=TRUE, then draw cats eyes}
  \item{alpha}{A degree of transparency (0=transparent ... 1 not transparent)}
  \item{adjust}{Allows smoothing of density histograms when plotting variables like height}
  
  \item{freq}{if TRUE, then plot frequencies (n * density)}
  \item{restrict}{Restrict the density to the observed max and min of the data}
  \item{xlim}{if not specified, will be .5 beyond the number of variables}
  \item{ylim}{If not specified, determined by the data}
  \item{add}{Allows overplotting}
  \item{col}{Allows for specification of colours.  The default for 2 groups is 
          blue and red, for more group levels, rainbows.}
 \item{pch}{The plot character for the mean is by default a small filled circle.  To not show the mean, use pch=NA }
 \item{scale}{If NULL, scale the widths by the square root of sample size, otherwise scale by the value supplied.}
 \item{legend}{If not NULL, draw a legend at c(topleft,topright,top,left,right) }
  \item{\dots}{Other graphic parameters}
}
\details{
Describe the data using a violin plot. Change alpha to modify the shading.  
The grp variable may be used to draw separate violin plots for each of multiple groups.

For relatively smallish data sets (< 500-1000), it is informative to also show the actual data points. This done with the dots=TRUE option. The jitter value is arbitrarily set to .05, but making it larger (say .1 or .2) will display more points.
}
\value{
The density (y axis) by value (x axis)  of the data  (for densityBy) or a violin plot for each variable (perhaps broken down by groups)
}
\author{William Revelle
}
\note{Formula input added July 12, 2020}

\seealso{  \code{\link{describe}}, \code{\link{describeBy}} and  \code{\link{statsBy}} for descriptive statistics and \code{\link{error.bars}} \code{\link{error.bars.by}} and  \code{\link{bi.bars}} and \code{\link{scatterHist}} for other graphic displays }
\examples{
violin(psychTools::bfi[4:8])
violin(SATV + SATQ ~ gender, data=sat.act, grp.name =cs(MV,FV,MQ,FQ)) #formula input
violinBy(psychTools::bfi,var=4:7,grp ="gender",grp.name=c("M","F"))
densityBy(SATV  ~ gender,data =sat.act,legend=1)  #formula input 
  
}
\keyword{ multivariate }
\keyword{ hplot }
