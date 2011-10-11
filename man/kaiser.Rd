\name{kaiser}
\alias{kaiser}

\title{Apply the Kaiser normalization when rotating factors}

\description{
Kaiser (1958) suggested normalizing factor loadings before rotating them, and then denormalizing them after rotation.  The GPArotation package does not (by default) normalize, nor does the \code{\link{fa}} function. Then, to make it more confusing, varimax in stats does,Varimax in GPArotation does not. \code{\link{kaiser}} will take the output of a non-normalized solution and report the normalized solution. 
}
\usage{
kaiser(f, rotate = "oblimin")
}
\arguments{
  \item{f}{A factor analysis output from \code{\link{fa}} or a factor loading matrix.
          }
  \item{rotate}{Any of the standard rotations avaialable in the GPArotation package. }
}
\details{Best results if called from an unrotated solution.  Repeated calls using a rotated solution will produce incorrect estimates of the correlations between the factors.
}

\value{See the values returned by GPArotation functions}
\references{
Kaiser, H. F. (1958) The varimax criterion for analytic rotation in factor analysis. Psychometrika 23, 187-200.
}
\author{
William Revelle
}
\note{
Prepared in response to a question about why \code{\link{fa}} oblimin results are different from SPSS.
}

\seealso{
\code{\link{fa}}
}
\examples{
f3 <- fa(Thurstone,3)
f3n <- kaiser(fa(Thurstone,3,rotate="none"))
factor.congruence(f3,f3n)
}

\keyword{ multivariate }
\keyword{ models}