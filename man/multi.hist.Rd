\name{multi.hist}
\alias{multi.hist}
\alias{histo.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Multiple histograms with density and normal fits on one page}
\description{Given a matrix or data.frame, produce histograms for each variable in a "matrix" form. Include normal fits and density distributions for each plot.

The number of rows and columns may be specified, or calculated.

May be used for single variables.
}
\usage{multi.hist(x,nrow=NULL, ncol=NULL,density=TRUE,main="Histogram, Density, and Normal Fit")
}
\arguments{
  \item{x}{ matrix or data.frame}
  \item{nrow}{number of rows in the plot}
  \item{ncol}{number of columns in the plot}
  \item{density}{density=TRUE, show the normal fits and density distributions}
  \item{main}{title for each panel}
}

\author{ William Revelle \cr
         Northwestern University \cr
         Evanston, Illinois \cr
 \email{  revelle@northwestern.edu           } \cr
 \url{http://personality-project.org/revelle.html}}

\examples{
#multi.hist(attitude[-1]) 
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
