\name{VSS.parallel}
\alias{VSS.parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare real and random VSS solutions}
\description{Another useful test for the number of factors is when the eigen values of a random matrix are greater than the eigen values of a a real matrix. Here we show VSS solutions to random data. 
}
\usage{
VSS.parallel(ncases, nvariables)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncases}{Number of simulated cases }
  \item{nvariables}{ number of simulated variables }
}

\value{VSS like output to be plotted by VSS.plot
}
\references{ ~put references to the literature/web site here ~ }
\author{ William Revelle}

\seealso{  \code{\link{VSS.plot}},  \code{\link{ICLUST}}, \code{\link{omega}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(ncases,nvariables)       #function call 
{  
 simdata=matrix(rnorm(ncases*nvariables),nrow=ncases,ncol=nvariables)    #make up simulated data
 testsim <- VSS(simdata,8,"none")
 return(testsim)
 }
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
