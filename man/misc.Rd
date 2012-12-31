\name{psych.misc}
\alias{psych.misc}
\alias{misc}
\alias{lowerCor}
\alias{lowerMat}
\alias{progressBar}
\alias{reflect}
\title{Miscellaneous helper functions for the psych package}
\description{This is a set of minor, if not trivial, helper functions.
lowerCor finds the correlation of x variables and then prints them using 
lowerMat which is a trivial, but useful, function to round off and print the lower triangle of a matrix.
reflect reflects the output of a factor analysis or principal components analysis so that one or more factors is reflected. (Requested by Alexander Weiss.)
progressBar prints out ...  as a calling routine (e.g., \code{\link{tetrachoric}} works through a tedious calculation.
}
\usage{
psych.misc()
lowerCor(x,digits=2,use="pairwise",method="pearson")
lowerMat(R, digits = 2)
reflect(f,flip=NULL)
progressBar(value,max,label=NULL) 
}
\arguments{
  \item{R}{A rectangular matrix or data frame (probably a correlation matrix)}
  \item{x}{A data matrix or data frame}
  \item{digits}{round to digits}
  \item{use}{Should pairwise deletion be done, or one of the other options to cor}
  \item{method}{"pearson", "kendall", "spearman"} 
  \item{value}{the current value of some looping variable}
  \item{max}{The maximum value the loop will achieve}
  \item{label}{what function is looping}
  \item{f}{The object returned from either a factor analysis (fa) or a principal components analysis (principal) }
  \item{flip}{The factor or components to be reversed keyed (by factor number)}
}
\value{ The lower triangle of a matrix, rounded to digits with titles abbreviated to digits + 3 (lowerMat) or a series of dots (progressBar).

lowerCor prints the lower diagonal correlation matrix but returns (invisibly) the full correlation matrix found with the use and method parameters.}

\examples{
lowerMat(Thurstone)
lowerCor(bfi[1:10])
f3 <- fa(Thurstone,3)
f3r <- reflect(f3,2)  #reflect the second factor
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}

