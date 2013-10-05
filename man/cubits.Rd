\name{cubits}
\alias{cubits}
\docType{data}
\title{Galton's example of the relationship between height and 'cubit' or forearm length}
\description{Francis Galton introduced the 'co-relation' in 1888 with a paper discussing how to measure the relationship between two variables.  His primary example was the relationship between height and forearm length.  The data table (cubits) is taken from Galton (1888).  Unfortunately, there seem to be some errors in the original data table in that the marginal totals do not match the table.

The data frame, \code{\link{heights}}, is converted from this table. 
}
\usage{data(cubits)}
\format{
  A data frame with 9 observations on the following 8 variables.
  \describe{
    \item{\code{16.5}}{Cubit length of lowest category}
    \item{\code{16.75}}{a numeric vector}
    \item{\code{17.25}}{a numeric vector}
    \item{\code{17.75}}{a numeric vector}
    \item{\code{18.25}}{a numeric vector}
    \item{\code{18.75}}{a numeric vector}
    \item{\code{19.25}}{a numeric vector}
    \item{\code{19.75}}{a numeric vector}
  }
}
\details{Sir Francis Galton (1888) published the first demonstration of the correlation coefficient.  The regression (or reversion to mediocrity) of the height to the length of the left forearm (a cubit) was found to .8.  There seem to be some errors in the table as published in that the row sums do not agree with the actual row sums. These data are used to create a matrix using \code{\link{table2matrix}} for demonstrations of analysis and displays of the data.
}
\seealso{  \code{\link{table2matrix}},  \code{\link{table2df}}, \code{\link{ellipses}}, \code{\link{heights}}, \code{\link{peas}},\code{\link{galton}}}
\source{Galton (1888)
}
\references{Galton, Francis (1888) Co-relations and their measurement. Proceedings of the Royal Society. London Series,45,135-145,
}
\examples{
data(cubits)
cubits
heights <- table2df(cubits,labs = c("height","cubit"))
ellipses(heights,n=1,main="Galton's co-relation data set")
ellipses(jitter(heights$cubit,3),jitter(heights$height,3),pch=".",
     main="Galton's co-relation data set") #add in some noise to see the points

}
\keyword{datasets}
