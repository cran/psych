\name{heights}
\alias{heights}
\docType{data}
\title{A data.frame of the Galton (1888) height and cubit data set.}
\description{Francis Galton introduced the 'co-relation' in 1888 with a paper discussing how to measure the relationship between two variables.  His primary example was the relationship between height and forearm length.  The data table (\code{\link{cubits}}) is taken from Galton (1888).  Unfortunately, there seem to be some errors in the original data table in that the marginal totals do not match the table.

The data frame, \code{\link{heights}}, is converted from this table using  \code{\link{table2df}}.
}
\usage{data(heights)}
\format{
  A data frame with 348 observations on the following 2 variables.
  \describe{
    \item{\code{height}}{Height in inches}
    \item{\code{cubit}}{Forearm length in inches}
  }
}
\details{Sir Francis Galton (1888) published the first demonstration of the correlation coefficient.  The regression (or reversion to mediocrity) of the height to the length of the left forearm (a cubit) was found to .8. The original table \code{\link{cubits}} is taken from Galton (1888). There seem to be some errors in the table as published in that the row sums do not agree with the actual row sums. These data are used to create a matrix using \code{\link{table2matrix}} for demonstrations of analysis and displays of the data.
}
\seealso{  \code{\link{table2matrix}},  \code{\link{table2df}}, \code{\link{cubits}}, \code{\link{ellipses}}, \code{\link{galton}}
}
\source{Galton (1888)
}
\references{Galton, Francis (1888) Co-relations and their measurement. Proceedings of the Royal Society. London Series,45,135-145,
}
\examples{
data(heights)
ellipses(heights,n=1,main="Galton's co-relation data set")

}
\keyword{datasets}
