\name{table2matrix}
\alias{table2matrix}
\alias{table2df}
\title{ Convert a table with counts to a matrix or data.frame representing those counts.}
\description{Some historical sets are reported as summary tables of counts in a limited number of bins.  Transforming these tables to data.frames representing the original values is useful for pedagogical purposes.  (E.g., transforming the original Galton table of height x cubits in order to demonstrate regression.) The column and row names must be able to be converted to numeric values.
}
\usage{
table2matrix(x, labs = NULL)
table2df(x, count=NULL,labs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A two dimensional table of counts with row and column names that can be converted to numeric values. }
  \item{count}{if present, then duplicate each row count times}
  \item{labs}{Labels for the rows and columns. These will be used for the names of the two columns of the resulting matrix }
}
\details{The original Galton (1888) of heights by cubits (arm length) is in tabular form. To show this as a correlation or as a scatter plot, it is useful to convert the table to a matrix or data frame of two columns.  

This function may also be used to convert an item response pattern table into a data table.  e.g., the Bock data set \code{\link{bock}}.  
}
\value{A matrix (or data.frame) of sum(x) rows and two columns.
}

\author{William Revelle}

\seealso{  \code{\link[psychTools]{cubits}} and  \code{\link{bock}} data sets}
\examples{
if(require(psychTools)) {
data(cubits)
cubit <- table2matrix(psychTools::cubits,labs=c("height","cubit"))
describe(cubit)
ellipses(cubit,n=1)
data(bock)
responses <- table2df(bock.table[,2:6],count=bock.table[,7],labs= paste("lsat6.",1:5,sep=""))
describe(responses)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

