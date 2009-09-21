\name{poly.mat}
\alias{poly.mat}
\title{Find polychoric correlations of item data}
\description{Uses John Fox's hetcor function (from polycor package) to find a matrix of polychoric correlations for integer data.  Essentially a wrapper for hetcor to convert integer item data into factor (categorical) data and then use hetcor.  Just a useful shortcut for subsequent factor analysis.
}
\usage{
poly.mat(x, short = TRUE, std.err = FALSE, ML = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix or data frame of integer data }
  \item{short}{ short=TRUE, just show the correlations, short=FALSE give the full hetcor output }
  \item{std.err}{std.err=FALSE does not report the standard errors (faster) }
  \item{ML}{ ML=FALSE  do a quick two step procedure, ML=TRUE, do longer maximum likelihood }
}
\details{Typical personality and item data are integer values (0,1 for ability;  1,2, 3, 4 for attitude scales).  The normal correlation procedures will find Pearson correlations (cor).  The polycor and hetcor functions from John Fox's polychor package will find polychoric correlations for categorical data.  This wrapper function converts integer data to categorical data and then calls hetcor.
}
\value{A matrix of polychoric correlations (if short=TRUE), otherwise a list of various estimates (see hetcor).
}
\author{ William Revelle}
\note{requires polycor}


\examples{
if(require(polycor)) {demo <- phi.demo()
round(demo$tetrachoric,2) #these should actually be all 1s but won't be
round(demo$phis,2)    #compare the phi (lower off diagonal and polychoric correlations (upper off diagonal)
} else {print("I am sorry, this demo requires the polycor package")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

