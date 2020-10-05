\name{scatterHist}
\alias{scatter.hist}
\alias{scatterHist}

\title{Draw a scatter plot with associated X and Y histograms, densities and correlation}
\description{Draw a X Y scatter plot with associated X and Y histograms with estimated densities.  Will also draw density plots by groups. Partly a demonstration of the use of layout. Also includes lowess smooth or linear model slope, as well as correlation.  Originally  from addicted to R example 78 with further modifications suggested by Jared Smith. 
}
\usage{scatterHist(x,y=NULL,smooth=TRUE,ab=FALSE,correl=TRUE,data=NULL,density=TRUE,ellipse=TRUE,  
  digits=2,method,cex.cor=1, title="Scatter plot + histograms", xlab=NULL,ylab=NULL,  
   smoother=FALSE,nrpoints=0, xlab.hist=NULL,ylab.hist=NULL,grid=FALSE, xlim=NULL, 
   ylim=NULL, x.breaks=11,y.breaks=11, x.space=0, y.space=0 ,freq=TRUE, x.axes=TRUE, 
   y.axes=TRUE,size=c(1,2),col=c("blue","red","black"),legend=NULL,alpha=.5,pch=21,...) 
 
 scatter.hist(x,y=NULL,smooth=TRUE,ab=FALSE,correl=TRUE,data=NULL,density=TRUE,
  ellipse=TRUE,  
  digits=2,method,cex.cor=1, title="Scatter plot + histograms", xlab=NULL,ylab=NULL,  
   smoother=FALSE,nrpoints=0, xlab.hist=NULL,ylab.hist=NULL,grid=FALSE, xlim=NULL, 
   ylim=NULL, x.breaks=11,y.breaks=11, x.space=0, y.space=0 ,freq=TRUE, x.axes=TRUE, 
   y.axes=TRUE,size=c(1,2),col=c("blue","red","black"),legend=NULL,alpha=.5,pch=21,...) 
}
\arguments{
  \item{x}{The X vector, or the first column of a  data.frame or matrix. Can be specified using formula input.  }
  \item{y}{The Y vector, of if X is a data.frame or matrix, the second column of X}
  \item{smooth}{if TRUE, then add a loess smooth to the plot}
  \item{ab}{if TRUE, then show the best fitting linear fit}
  \item{correl}{TRUE: Show the correlation}
  \item{data}{if using formula input, the data must be specified}
   \item{density}{TRUE: Show the estimated densities}
  \item{ellipse}{TRUE: draw 1 and 2 sigma ellipses and smooth}
  \item{digits}{How many digits to use if showing the correlation}
  \item{method}{Which method to use for correlation ("pearson","spearman","kendall") 
    defaults to "pearson"}
  \item{smoother}{if TRUE, use smoothScatter instead of plot. Nice for large N.}
  \item{nrpoints}{If using smoothScatter, show nrpoints as dots. Defaults to 0}
 \item{grid}{If TRUE, show a grid for the scatter plot.}
  \item{cex.cor}{Adjustment for the size of the correlation}
  \item{xlab}{Label for the x axis}
  \item{ylab}{Label for the y axis}
  \item{xlim}{Allow specification for limits of x axis, although this seems to just work for the scatter plots.}
  \item{ylim}{Allow specification for limits of y axis} 
  \item{x.breaks}{Number of breaks to suggest to the x axis histogram.}
  \item{y.breaks}{Number of breaks to suggest to the y axis histogram.}
  \item{x.space}{space between bars}
  \item{y.space}{Space between y bars}
  \item{freq}{Show frequency counts, otherwise show density counts}
  \item{x.axes}{Show the x axis for the x histogram}
  \item{y.axes}{Show the y axis for the y histogram}
  \item{size}{The sizes of the ellipses (in sd units).  Defaults to 1,2}
  \item{col}{Colors to use when showing groups}
  \item{alpha}{Amount of transparency in the density plots}
  \item{legend}{Where to put a legend  c("topleft","topright","top","left","right")}
  \item{pch}{Base plot character (each group is one more)}
  \item{xlab.hist}{Not currently available}
  \item{ylab.hist}{Label for y axis histogram.  Not currently available.} 
  \item{title}{An optional title}
  \item{\dots}{Other parameters for graphics}
}
\details{Just a straightforward application of layout and barplot, with some tricks taken from \code{\link{pairs.panels}}.  The various options allow for correlation ellipses (1 and 2 sigma from the mean), lowess smooths, linear fits, density curves on the histograms, and the value of the correlation.  ellipse = TRUE implies smooth = TRUE.  The grid option provides a background grid to the scatterplot.

Formula input allows specification of grouping variables as well.  )
}
\author{William Revelle}
\note{
Adapted from Addicted to R example 78.  Modified following some nice suggestions from Jared Smith. 
}

\seealso{\code{\link{pairs.panels}} for multiple plots, \code{\link{multi.hist}} for multiple histograms.
}
\examples{
data(sat.act)
with(sat.act,scatterHist(SATV,SATQ))
scatterHist(SATV ~ SATQ,data=sat.act)  #formula input

#or for something a bit more splashy
scatter.hist(sat.act[5:6],pch=(19+sat.act$gender),col=c("blue","red")[sat.act$gender],grid=TRUE)
#better yet
scatterHist(SATV ~ SATQ + gender,data=sat.act) #formula input with a grouping variable
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ hplot }
