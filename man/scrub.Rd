\name{scrub}
\alias{scrub}

\title{A utility for basic data cleaning and recoding.  Changes values outside of minimum and maximum limits to NA.}
\description{A tedious part of data analysis is addressing the problem of miscoded data that needs to be converted to NA.  For a given data.frame or matrix, scrub will set all values of columns from=from to to=to that are less than a set (vector) of min values or more than a vector of max values to NA. Can also be used to do basic recoding of data for all values=isvalue to newvalue.
}
\usage{
scrub(x, where, min, max,isvalue,newvalue)
}

\arguments{
  \item{x}{a data frame or matrix}
  \item{where}{The variables to examine.  (Can be by name or by column number)}
  \item{min}{a vector of minimum values that are acceptable}
  \item{max}{a vector of maximum values that are acceptable}
 \item{isvalue}{a vector of values to be converted to newvalue (one per variable)}
\item{newvalue}{a vector of values to replace those that match isvalue}}
\details{Solves a tedious problem that can be done directly but that is sometimes awkward.  Will either replace specified values with NA or 
}
\value{
The corrected data frame.}

\author{William Revelle
}
\note{Probably could be optimized to avoid one loop
}

\seealso{ \code{\link{reverse.code}}, \code{\link{rescale}} for other simple utilities.
}
\examples{
data(attitude)
x <- scrub(attitude,isvalue =c(30,40,50), newvalue = c(930,940,950))  #will change all occurences
x1 <- scrub(attitude, where=4, isvalue =c(30,40,50), newvalue = c(930,940,950)) #will just do it for the 4th column
#get rid of a complicated set of cases and replace with missing values
y <- scrub(attitude,where=2:4,min=c(20,30,40),max= c(120,110,100),isvalue= c(32,43,54))
y
scrub(attitude,where="learning",isvalue=68,newvalue=999) #change a column by name
scrub(attitude,where="learning",min=45,newvalue=999) #change a column by name

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
%code is in reverse.code
