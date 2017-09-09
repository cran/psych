\name{residuals.psych}
\alias{residuals.psych}
\alias{resid.psych}
\title{Extract residuals from various psych objects}

\description{Residuals in the various psych functions are extracted and then may be "pretty" printed.}

\usage{
\method{residuals}{psych}(object,diag=TRUE,...)
\method{resid}{psych}(object,diag=TRUE,...)
}

\arguments{
  \item{object}{The object returned by a psych function.}
  \item{diag}{if FALSE, then convert the diagonal of the residuals to NA}
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
sum(residuals(f3)^2) #include diagonal
sum(residuals(f3,diag=FALSE)^2,na.rm=TRUE) #drop diagonal
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ models}