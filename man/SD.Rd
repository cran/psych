\name{SD}
\alias{SD}
\title{ Find the Standard deviation for a vector, matrix, or data.frame - do not return error if there are no cases }
\description{Find the standard deviation of a vector, matrix, or data.frame.  In the latter two cases, return the sd of each column.  Unlike the sd function, return NA if there are no observations rather than throw an error.
}
\usage{
SD(x, na.rm = TRUE)   #deprecated
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, data.frame, or matrix }
  \item{na.rm}{na.rm is assumed to be TRUE}
}
\details{Finds the standard deviation of a vector, matrix, or data.frame.  Returns NA if no cases.

Just an adaptation of the stats:sd function to return the functionality found in R < 2.7.0 or R >= 2.8.0
Because this problem seems to have been fixed, SD will be removed eventually.}
\value{The standard deviation
}

\author{ William Revelle }
\note{ Until R 2.7.0, sd would return a NA rather than an error if no cases were observed.  SD brings back that functionality.  Although unusual, this condition will arise when analyzing data with high rates of missing values.
This function will probably be removed as 2.7.0 becomes outdated.}

\seealso{ These functions use SD rather than sd:  \code{\link{describe.by}}, \code{\link{skew}}, \code{\link{kurtosi}}  }

\examples{
data(attitude)
sd(attitude) #all complete
attitude[,1] <- NA
SD(attitude) #missing a column
describe(attitude)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

