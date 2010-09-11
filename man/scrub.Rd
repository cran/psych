\name{scrub}
\alias{scrub}

\title{A utility for basic data cleaning.  Changes values outside of minimum and maximum limits to NA
}
\description{A tedious part of data analysis is addressing the problem of miscoded data that needs to be converted to NA.  For a given data.frame or matrix, scrub will set all values of columns from=from to to=to that are less than a set (vector) of min values or more than a vector of max values to NA.
}
\usage{
scrub(x, where, min, max,isvalue)
}

\arguments{
  \item{x}{a data frame or matrix
}
  \item{where}{The variables to examine.  (Can be by name or by column number)
}

  \item{min}{
a vector of minimum values that are acceptable
}
  \item{max}{
a vector of maximum values that are acceptable
}
 \item{isvalue}{
a vector of values to be converted to NA (one per variable)
}
}
\details{
Solves a tedious problem that can be done directly but that is sometimes awkward.
}
\value{
The corrected data frame.
}

\author{William Revelle
}
\note{Probably could be optimized to avoid one loop
}

\seealso{ \code{\link{reverse.code}}, \code{\link{rescale}} for other simple utilities.
}
\examples{
x <- matrix(1:150,ncol=10,byrow=TRUE)
#get rid of a complicated set of cases
y <- scrub(x,where=2:4,min=c(20,30,40),max= c(120,110,100),isvalue= c(32,43,54))
y
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
%code is in reverse.code
