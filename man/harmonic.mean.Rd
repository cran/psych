\name{harmonic.mean}
\alias{harmonic.mean}
\title{ Find the harmonic mean of a vector, matrix, or columns of a data.frame}
\description{The harmonic mean is merely the reciprocal of the arithmetic mean of the reciprocals.
}
\usage{
harmonic.mean(x,na.rm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector, matrix, or data.frame }
  \item{na.rm}{na.rm=TRUE remove NA values before processing}
}
\details{Included as an example for teaching about functions. As well as for a discussion of how to estimate central tendencies.
}
\value{
  The harmonic mean(s)
}

\note{Included as a simple demonstration of how to write a function}

\examples{
x <- seq(1,5)
x2 <- x^2
x2[2] <- NA
X <- data.frame(x,x2)
harmonic.mean(x)
harmonic.mean(x2)
harmonic.mean(X)
harmonic.mean(X,FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
