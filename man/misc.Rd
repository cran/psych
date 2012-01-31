\name{psych.misc}
\alias{psych.misc}
\alias{misc}
\alias{lower.mat}
\alias{progressBar}
\title{Miscellaneous helper functions for the psych package}
\description{This is a set of minor, if not trivial, helper functions.

lower.mat is a trivial, but useful,  function to round off and print the lower triangle of a matrix.

progressBar prints out ...  as it works through a tedious calculation.
}
\usage{
psych.misc()
lower.mat(R, digits = 2)
progressBar(value,max,label=NULL) 
}
\arguments{
  \item{R}{A matrix or data frame}
  \item{digits}{round to digits}
  \item{value}{the current value of some looping variable}
  \item{max}{The maximum value the loop will achieve}
  \item{label}{what function is looping}
}
\value{ The lower triangle of a matrix, rounded to digits with titles abbreviated to digits + 3 (lower.mat) or a series of dots.}

\examples{
lower.mat(Thurstone)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

