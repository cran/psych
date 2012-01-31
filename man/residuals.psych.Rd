\name{residuals.psych}
\alias{residuals.psych}
\alias{resid.psych}
\title{Extract residuals from various psych objects}

\description{Residuals in the various psych functions are extracted and then may be "pretty" printed.}

\usage{
\method{residuals}{psych}(object,...)
\method{resid}{psych}(object,...)
}

\arguments{
  \item{object}{The object returned by a psych function.}
  \item{...}{Other parameters to be passed to residual (ignored but required by the generic function)}
}
\details{
Currently implemented for \code{\link{fa}}, \code{\link{principal}},   \code{\link{omega}}, \code{\link{irt.fa}}, and \code{\link{fa.extension}}. 
}
\value{
residuals:  a matrix of residual estimates}

\author{William Revelle
}
\examples{
f3 <- fa(Thurstone,3)
residuals(f3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}