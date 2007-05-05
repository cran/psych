\name{matrix.addition}
\alias{\%+\%}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to add two vectors or matrices }
\description{It is sometimes convenient to add two vectors or matrices in an operation analogous to matrix multiplication. For matrices nXm and mYp, the matrix sum  of the i,jth element of nSp = sum(over m) of iXm + mYj. 
}
\usage{
x \%+\% y }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a n by m matrix (or vector if m= 1)}
  \item{y}{ a m by p matrix (or vector if m = 1)}
}
\details{Used in such problems as Thurstonian scaling.  Although not technically matrix addition, as pointed out by Krus, there are many applications where the sum or difference of two vectors or matrices is a useful operation.  This can be done, of course, through 
}
\value{a n by p matix of sums

}
\references{Krus, D. J. (2001) Matrix addition. Journal of Visual Statistics, 1, (February, 2001).}
\author{William Revelle}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
x <- seq(1,4)
z <- x \%+\% -t(x)
x
z
x <- matrix(seq(1,6),ncol=2)
y <- matrix(seq(1,10),nrow=2)
z <- x \%+\% y
x
y
z
## The function is currently defined as
"\%+\%" <- function(x,y) { 
 if(!is.matrix(x)) {
	if(is.vector(x)) {x <- as.matrix(x)} else stop("x must be either a vector or a matrix")}
if(!is.matrix(y)) {
	if(is.vector(y)) {y <- as.matrix(y)} else stop("y must be either a vector or a matrix")}
n.x <- dim(x)[1]
n.y <- dim(y)[2]
n.k <- dim(x)[2]
if (n.k != dim(y)[1]) {warning("Matrices should be comparable")}
#first find sum vectors
x <- rowSums(x)
y <- colSums(y)
one <- as.vector(rep(1,n.y))  #to duplicate x n.y times
one.y <- as.vector(rep(1,n.x)) #to duplicate y n.x times
xy  <- x %*% t(one) + t(y %*% t(one.y) )  #sum the vectors in a rectangular array
  return(xy) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
