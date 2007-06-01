\name{harmonic.mean}
\alias{harmonic.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the harmonic mean of a vector, matrix, or columns of a data.frame}
\description{The harmonic mean is merely the reciprocal of the arithmetic mean of the reciprocals.
}
\usage{
harmonic.mean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector, matrix, or data.frame }
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
harmonic.mean(x)
harmonic.mean(x2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
