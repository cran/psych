\name{skew}
\alias{skew}
\alias{kurtosi}

\title{ Calculate skew or kurtosis for a vector, matrix, or data.frame}
\description{Find the skew and kurtosis for each variable in a data.frame or matrix.  Unlike skew and kurtosis in e1071, this calculates a different skew for each variable or column of a data.frame/matrix.
}
\usage{
skew(x, na.rm = TRUE)
kurtosi(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data.frame or matrix }
  \item{na.rm}{ how to treat missing data }
}
\details{given a matrix or data.frame x, find the skew or kurtosis for each column.
}
\value{if input is a matrix or data.frame, skew (kurtosi) is a vector of skews (kurtosi)}
\author{William Revelle }

\note{The mean function supplies means for the columns of a data.frame, but the overall mean for a matrix.  Mean will throw a warning for non-numeric data, but colMeans stops with non-numeric data. Thus, the function uses either mean (for data frames) or colMeans (for matrices).  This is true for skew and kurtosi as well.

}

\seealso{\code{\link{describe}}, \code{\link{describe.by}},   }


\examples{
round(skew(attitude),2)
round(kurtosi(attitude),2)

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
