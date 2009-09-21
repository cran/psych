\name{tr}
\alias{tr}

\title{Find the trace of a square matrix}
\description{Hardly worth coding, if it didn't appear in so many formulae in psychometrics, the trace of a (square) matrix is just the sum of the diagonal elements. 
}
\usage{
tr(m)
}

\arguments{
  \item{m}{A square matrix }
}
\details{The tr function is used in various matrix operations and is the sum of the diagonal elements of a matrix.
}
\value{ The sum of the diagonal elements of a square matrix.  \cr
i.e. tr(m)  <- sum(diag(m)). 
}




\examples{
 m <- matrix(1:16,ncol=4)
 m
 tr(m)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

