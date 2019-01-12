\name{cattell}
\alias{cattell}
\docType{data}
\title{12 cognitive variables from Cattell (1963) }
\description{
Rindskopf and Rose (1988) use this data set to demonstrate confirmatory second order factor models.  It is a nice example data set to explore hierarchical structure and alternative factor solutions. It contains measures of fluid and crystallized intelligence.
}
\usage{data("cattell")}
\format{
  A correlation matrix of the following 12 variables from 277 7th and 8th graders
  \describe{
\item{Verbal}{A verbal ability test from Thurstone}
\item{Verbal2}{A verbal ability test from Thurstone}
\item{Space1}{A Spatial ability test from Thurstone}
\item{Space2}{A Spatial ability test from Thurstone}
\item{Reason1}{A reasoning test from Thurstone}
\item{Reason2}{A reasoning  test from Thurstone}
\item{Number1}{A Numerical ability test from Thurstone}
\item{Number2}{A Numerical ability test from Thurstone}
\item{IPATSer}{A "culture fair"  series from the IPAT}
\item{IPATCLAS}{A "culture fair" classification test from the IPAT}
\item{IPATMatr}{A "culture fair" matrix reasoning  test from the IPAT}
\item{IPATTop}{A "culture fair" topology test from the IPAT}
}
}
\details{
Cattell (1963) reported on 8 cognitive variables from Thurstone and four from the Institute for Personality Assessment Test (IPAT).  Rindskopf and Rose (1988) use this data set as an example of second order factor analysis. It is thus a nice set for examining alternative solutions such as bifactor rotation, \code{\link{omega}} hierarchical, as well as \code{\link{esem}} and \code{\link{interbattery}} factor analysis.
}
\source{
David Rindskopf and Tedd Rose, (1988) Some Theory and Applications of Confirmatory Second- Order Factor Analysis, Multivariate Behavioral Research, 23, 51-67.}
\references{

Cattell, R. B. (1963).Theory of fluid and crystallized intelligence: A critical experiment. Journal of Educational Psychology, 54, 1-22.

David Rindskopf and Tedd Rose, (1988) Some Theory and Applications of Confirmatory Second- Order Factor Analysis, Multivariate Behavioral Research, 23, 51-67.
}
\examples{
data(cattell)
corPlot(cattell,numbers=TRUE,upper=FALSE,diag=FALSE,
             main="12 cognitive variables from Cattell (1963)",xlas=2)
}
\keyword{datasets}
