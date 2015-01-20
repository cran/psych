\name{Harman.5}
\alias{Harman.5}
\docType{data}
\title{5 socio-economic variables from Harman (1967)
}
\description{
Harman (1967) uses 5 socio-economic variables for demonstrations of principal components and factor analysis.  This example is used in the SAS manual for Proc Factor as well.}
\usage{data(Harman.5)}
\format{
  A data frame with 12 observations on the following 5 variables.
  \describe{
    \item{\code{population}}{a numeric vector}
    \item{\code{schooling}}{a numeric vector}
    \item{\code{employment}}{a numeric vector}
    \item{\code{professional}}{a numeric vector}
    \item{\code{housevalue}}{a numeric vector}
  }
}
\details{
Harman reports that the data "were taken (not entirely arbitrarily) from a study of the Los Angeles Standard Metropolitan Statistical Area.  The twelve individuals are used in the examples are census tracts." (p 13).
}
\source{
Harman, Harry Horace (1967), Modern factor analysis. Chicago, University of Chicago Press. }
\references{
SAS users manual, chapter 26: pages 1123-1192}
\examples{
data(Harman.5)
if(require('GPArotation')){
pc2 <- principal(Harman.5,2,scores=TRUE)
pc2$residual
biplot(pc2,main="Biplot of the Harman 5 socio-demographic variables") }
}
\keyword{datasets}
