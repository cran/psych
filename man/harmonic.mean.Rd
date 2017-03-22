\name{harmonic.mean}
\alias{harmonic.mean}
\title{ Find the harmonic mean of a vector, matrix, or columns of a data.frame}
\description{The harmonic mean is merely the reciprocal of the arithmetic mean of the reciprocals.
}
\usage{
harmonic.mean(x,na.rm=TRUE,zero=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector, matrix, or data.frame }
  \item{na.rm}{na.rm=TRUE remove NA values before processing}
  \item{zero}{If TRUE, then if there are any zeros, return 0, else, return the harmonic mean of the non-zero elements}
}
\details{Included as an example for teaching about functions. As well as for a discussion of how to estimate central tendencies.
Also used in \code{\link{statsBy}} to weight by the harmonic mean.

Values of 0 can be included (in which case the harmonic.mean = 0) or converted to NA according to the zero option.  

Added the zero option, March, 2017.
}
\value{
  The harmonic mean(s)
}

\note{Included as a simple demonstration of how to write a function}

\examples{
x <- seq(1,5)
x2 <- x^2
x2[2] <- NA
y <- x - 1
X <- data.frame(x,x2,y)
harmonic.mean(x)
harmonic.mean(x2)
harmonic.mean(X)
harmonic.mean(X,na.rm=FALSE)
harmonic.mean(X,zero=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
