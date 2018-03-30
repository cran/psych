\name{dfOrder}
\alias{dfOrder}

\title{Sort (order) a dataframe or matrix by multiple columns
}
\description{Although \code{\link{order}} will order a vector, and it is possible to order several columns of a data.frame by specifying each column individually in the call to order, \code{\link{dfOrder}} will order a dataframe by as many columns as desired.
}
\usage{
dfOrder(object, columns,absolute=FALSE)
}

\arguments{
  \item{object}{The data.frame to be sorted}

  \item{columns}{Column numbers to use for sorting.  If positive, then they will be sorted in increasing order. If negative, then in decreasing order}
  \item{absolute}{If TRUE, then sort the absolute values}
}
\details{
This is just a simple helper function to reorder data.frames.  Originally developed to organize IRT  output from the ltm package. It is a basic add on to the order function.  

(Completely rewritten for version 1.8.1.) 
}
\value{
The original data frame is now in sorted order.
}

\author{William Revelle
}


\seealso{ Other useful file manipulation functions include \code{\link{read.file}} to read in data from a file or \code{\link{read.clipboard}} from the clipboard,  \code{\link{fileScan}}, \code{\link{filesList}}, \code{\link{filesInfo}},  and \code{\link{fileCreate}}

\code{\link{dfOrder}} code is used in the  \code{\link{test.irt}} function to combine ltm and \code{\link{sim.irt}} output.

}
\examples{
x <- data.frame(matrix(sample(1:4,64,replace=TRUE),ncol=4))
dfOrder(x)  # sort by all columns
dfOrder(x,c(1,4))  #sort by the first and 4th column
dfOrder(x,c(1,-2))  #sort by the first in increasing order, 
   #the second in decreasing order

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip }% use one of  RShowDoc("KEYWORDS")
\keyword{utilities }% __ONLY ONE__ keyword per line
