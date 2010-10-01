\name{plot.psych}
\alias{plot.psych}
\alias{plot.poly}
\alias{plot.irt}




\title{Plotting functions for the psych package of class ``psych"}
\description{Combines several plotting functions into one for objects of class ``psych".  This can be used to plot the results of \code{\link{fa}}, \code{\link{irt.fa}}, \code{\link{VSS}}, \code{\link{ICLUST}}, \code{\link{omega}}, \code{\link{factor.pa}}, or \code{\link{principal}}. 
}
\usage{
plot.psych(x,labels=NULL,...)
plot.irt(x,xlab,ylab,main,D,type=c("ICC","IIC","test"),...)
plot.poly(x,D,xlab,ylab,ylim,main,type=c("ICC","IIC","test"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The object to plot }
  \item{labels}{Variable labels}
  \item{xlab}{Label for the x axis  -- defaults to Latent Trait}
  \item{ylab}{Label for the y axis}
  \item{ylim}{Specify the limits for the y axis}
  \item{main}{Main title for graph}
  \item{type}{ICC plots items, IIC plots item information, test plots test information, defaults to IIC for plot.irt, to ICC for plot.poly}
  \item{D}{The discrimination parameter}
  \item{...}{other calls to plot}
}
\details{Passes the appropriate values to plot.  For plotting the results of \code{\link{irt.fa}}, there are three options: type = "ICC" (default) will plot the item characteristic respone function.  type = "IIC" will plot the item information function, and type= "test" will plot the test information function.

These are calls to the generic plot function that are intercepted for objects of type "psych".  More precise plotting control is available in the separate plot functions.  plot may be used for psych objects returned from \code{\link{fa}}, \code{\link{irt.fa}}, \code{\link{ICLUST}}, \code{\link{omega}}, as well as  \code{\link{principal}} 

}
\value{Graphic output for factor analysis, cluster analysis and item response analysis. 

 
}

\author{William Revelle }
\note{ More precise plotting control is available in the separate plot functions.

 
}
\seealso{ \code{\link{VSS.plot}} and \code{\link{factor.plot}}, \code{\link{cluster.plot}}, \code{\link{fa}}, \code{\link{irt.fa}}, \code{\link{VSS}}, \code{\link{ICLUST}}, \code{\link{omega}}, \code{\link{factor.pa}}, or \code{\link{principal}}   
}
\examples{
test.data <- Harman74.cor$cov
f4 <- fa(test.data,4)
plot(f4)

#not run
#data(bfi)
#e.irt <- irt.fa(bfi[11:15])  #just the extraversion items
#plot(e.irt)   #the information curves
#
#ic <- iclust(test.data,3)   #shows hierarchical structure 
#plot(ic)                    #plots loadings
#


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

