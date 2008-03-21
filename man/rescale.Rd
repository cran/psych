\name{rescale}
\alias{rescale}
\title{Function to convert scores to ``conventional
" metrics}
\description{Psychologists frequently report data in terms of transformed scales such as ``IQ" (mean=100, sd=15, ``SAT/GRE" (mean=500, sd=100), ``ACT" (mean=18, sd=6), ``T-scores" (mean=50, sd=10), or ``Stanines" (mean=5, sd=2). The \code{\link{rescale}} function converts the data to standard scores and then rescales to the specified mean and standard deviation}
\usage{
rescale(x, mean = 100, sd = 15,df=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix or data frame }
  \item{mean}{Desired mean of the rescaled scores}
  \item{sd}{Desired standard deviation of the rescaled scores}
  \item{df}{if TRUE, returns a data frame, otherwise a matrix}
}


\value{A data.frame or matrix of rescaled scores.
  
}

\author{ William Revelle }

\seealso{  See Also  \code{\link{scale}} }
\examples{
T <- rescale(attitude,50,10)
describe(T)
}

\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{univar}
