\name{dfOrder}
\alias{dfOrder}

\title{Sort (order) a dataframe or matrix by multiple columns
}
\description{Although order will order a vector, and it is possible to order several columns of a data.frame by specifying each column individually in the call to order, dfOrder will order a dataframe by as many columns as desired.
}
\usage{
dfOrder(object, columns)
}

\arguments{
  \item{object}{The data.frame to be sorted}

  \item{columns}{Column numbers to use for sorting.  If positive, then they will be sorted in increasing order. If negative, then in decreasing order}
}
\details{
This is just a simple helper function to reorder data.frames.  Developed to organize IRT  output from the ltm package. It is a basic add on to the order function. 
}
\value{
The original data frame is now in sorted order.
}

\author{William Revelle
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
This code is used in the  \code{\link{test.irt}} function to combine ltm and \code{\link{sim.irt}} output.
}
\examples{
x <- data.frame(matrix(sample(1:4,50,replace=TRUE),ncol=5))
dfOrder(x,c(1,4))
dfOrder(x,c(1,-2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip }% use one of  RShowDoc("KEYWORDS")
\keyword{utilities }% __ONLY ONE__ keyword per line
