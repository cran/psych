\name{densityBy}
\alias{densityBy}
\alias{violinBy}

\title{Create a 'violin plot' or density plot of the distribution of a set of variables}

\description{Among the many ways to describe a data set, one is density plot or violin plot of the data.  This is similar to a box plot but shows the actual distribution.
Median and 25th and 75th percentile lines are added to the display. If a grouping variable is specified, densityBy will draw violin plots for each variable and for each group. }
\usage{
densityBy(x,grp=NULL,grp.name=NULL,ylab="Observed",xlab="",main="Density plot",density=20,
         restrict=TRUE,xlim=NULL,add=FALSE,col=NULL,pch=20, ...)
violinBy(x,grp=NULL,grp.name=NULL,ylab="Observed",xlab="",main="Density plot",density=20,
         restrict=TRUE,xlim=NULL,add=FALSE,col=NULL,pch=20, ...) 
}


\arguments{
  \item{x}{A matrix or data.frame}
  \item{grp}{A grouping variable}
  \item{grp.name}{If the grouping variable is specified, the what names should be give to the group? Defaults to 1:ngrp}
  \item{ylab}{The y label}
  \item{xlab}{The x label}
  \item{main}{Figure title}
  \item{density}{How many lines per inch to draw}
  \item{restrict}{Restrict the density to the observed max and min of the data}
  \item{xlim}{if not specified, will be .5 beyond the number of variables}
  \item{add}{Allows overplotting}
  \item{col}{Allows for specification of colours.  The default for 2 groups is 
          blue and red, for more group levels, rainbows.}
 \item{pch}{The plot character for the mean is by default a small filled circle.  To not show the mean, use pch=NA }
  \item{\dots}{Other graphic parameters}
}
\details{
Describe the data using a violin plot. Change density to modify the shading.  density=NULL will fill with col.
The grp variable may be used to draw separate violin plots for each of multiple groups.
}
\value{
The density plot of the data.
}
\author{William Revelle
}
\note{Nothing yet}
\seealso{  \code{\link{describe}}, \code{\link{describeBy}} and  \code{\link{statsBy}} for descriptive statistics and \code{\link{error.bars}} and \code{\link{error.bars.by}} for graphic displays }
\examples{
densityBy(bfi[1:5])
#not run
#violinBy(bfi[1:5],grp=bfi$gender,grp.name=c("M","F"))
#densityBy(sat.act[5:6],sat.act$education,col=rainbow(6))
  
}
\keyword{ multivariate }
\keyword{ hplot }
