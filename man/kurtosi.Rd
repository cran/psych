\name{kurtosi}
\alias{kurtosi}

\title{ Kurtosis of a vector, matrix, or data frame}
\description{Find the kurtosis of a vector, matrix, or dataframe.
}
\usage{
kurtosi(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector, matrix, or data frame }
  \item{na.rm}{ na.rm =TRUE removes missing data from the column  }
}
\details{Kurtosis in the E1071 package finds the kurtosis for a single vector.  This does it for matrices and dataframes. Used in the describe function. 


}
\value{
  \item{kurtosi}{a vector of the kurtosis for each column of the matrix}
 
}


\author{William Revelle }

\note{The mean function supplies means for the columns of a data.frame, but the overall mean for a matrix.  Mean will throw a warning for non-numeric data, but colMeans stops with non-numeric data. Thus, the function uses either mean (for data frames) or colMeans (for matrices).  This is true for skew and kurtosi as well.}
\seealso{  \code{\link{skew}}, \code{\link{describe}} }
\examples{

 round(kurtosi(attitude),2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

