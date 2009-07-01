\name{VSS.parallel}
\alias{VSS.parallel}
\title{Compare real and random VSS solutions}
\description{Another useful test for the number of factors is when the eigen values of a random matrix are greater than the eigen values of a a real matrix. Here we show VSS solutions to random data. 
}
\usage{
VSS.parallel(ncases, nvariables,scree=FALSE,rotate="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncases}{Number of simulated cases }
  \item{nvariables}{ number of simulated variables }
  \item{scree}{Show a scree plot for random data -- see  \code{\link{omega}}}
  \item{rotate}{rotate="none" or rotate="varimax"}
}

\value{VSS like output to be plotted by VSS.plot
}
\references{Very Simple Structure (VSS)}
\author{ William Revelle}

\seealso{  \code{\link{fa.parallel}}, \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}} }
\examples{
#VSS.plot(VSS.parallel(200,24))
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
