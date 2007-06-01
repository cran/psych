\name{VSS.simulate}
\alias{VSS.simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ create VSS like data}
\description{Simulation is one of most useful techniques in statistics and psychometrics.  Here we simulate a correlation matrix with a simple structure composed of a specified number of factors.  Each item is assumed to have complexity one.
}
\usage{
VSS.simulate(ncases=1000, nvariables=16, nfactors=4, meanloading=.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncases}{ number of simulated subjects }
  \item{nvariables}{ Number of variables }
  \item{nfactors}{ Number of factors to generate }
  \item{meanloading}{with a mean loading }
}

\value{a ncases x nvariables  matrix 
}

\author{ William Revelle}

\seealso{ \code{\link{VSS}}, \code{\link{ICLUST}} }
\examples{
\dontrun{
simulated <- VSS.simulate(1000,20,4,.6)
vss <- VSS(simulated,rotate="varimax")
VSS.plot(vss)
}

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
