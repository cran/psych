\name{plot.psych}
\alias{plot.psych}

\title{Plotting functions for the psych package of class ``psych"}
\description{Combines several plotting functions into one for objects of class ``psych".  This can be used to plot the results of \code{\link{VSS}}, \code{\link{ICLUST}}, \code{\link{omega}}, \code{\link{factor.pa}}, or \code{\link{principal}}. 
}
\usage{
plot.psych(x,labels=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The object to plot }
  \item{labels}{Variable labels}
  \item{...}{other calls to plot}
}
\details{Passes the appropriate values to plot
}
\value{
 
}

\author{William Revelle }
\note{ More precise plotting control is available in the separate plot functions.

 
}
\seealso{ \code{\link{VSS.plot}} and \code{\link{cluster.plot}} }
\examples{
test.data <- Harman74.cor$cov
f4 <- factor.pa(test.data,4)
plot(f4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

