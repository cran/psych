\name{describe}
\alias{describe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Basic descriptive statistics useful for psychometrics }
\description{
There are many summary statistics available in R; this function
  provides the ones most useful for scale construction and item analysis in classic psychometrics. 
  Range is most useful for the first pass in a data set, to check for coding errors. 
}
\usage{
describe(x, digits = 2, na.rm = TRUE, skew = TRUE, ranges = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame or matrix}
  \item{digits}{ How many significant digits to report}
  \item{na.rm}{The default is to delete missing data }
  \item{skew}{ Should the skew and kurtosis be calculated? }
  \item{ranges}{ Should the range be calculated? }
}
\details{In basic data analysis it is vital to get basic descriptive statistics. 
Procedures such as \code{\link{summary}} and hmisc::describe do so.  The describe function in the \code{\link{psych}} package is meant to produce the most frequently requested stats in psychometric and psychology studies, and to produce them in an easy to read data.frame. The results from describe can be used in graphics functions (e.g., \code{\link{error.crosses}}).

The range statistics (min, max, range) are most useful for data checking to detect coding errors, and should be found in early analyses of the data.  

Although describe will work on data frames as well as matrices, it is important to realize that for data frames, descriptive statistics will be reported only for those variables where this makes sense (i.e., not for factors or for alphanumeric data). 

In a typical study, one might read the data in from the clipboard (\code{\link{read.clipboard}}), show the splom plot of the correlations (\code{\link{pairs.panels}}), and then describe the data. 


}
\value{
 A data.frame of the relevant statistics: \cr
   item name \cr
   item number \cr
   number of valid cases\cr
   mean\cr
   standard deviation\cr
   median\cr
   mad: median absolute deviation (from the median) \cr
   minimum\cr
   maximum\cr
   skew\cr
   kurtosis\cr
   standard error\cr
}

\note{Describe uses either the mean or colMeans functions depending upon whether the data are a data.frame or a matrix. The mean function supplies means for the columns of a data.frame, but the overall mean for a matrix.  Mean will throw a warning for non-numeric data, but colMeans stops with non-numeric data. Thus, the describe function uses either mean (for data frames) or colMeans (for matrices).  This is true for skew and kurtosi as well.}
\author{ 
 
 \url{http://personality-project.org/revelle.html} \cr

Maintainer: William Revelle \email{revelle@northwestern.edu} \cr
}
\seealso{ \code{\link{describe.by}}, \code{\link{skew}}, \code{\link{kurtosi}}, \code{\link{pairs.panels}}, \code{\link{read.clipboard}}, \code{\link{error.crosses}} }
\examples{

describe(attitude)

describe(attitude,skew=FALSE)   #attitude is taken from R data sets






}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{univar}
