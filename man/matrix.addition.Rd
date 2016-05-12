\name{matrix.addition}
\alias{matrix.addition}
\alias{\%+\%}
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
\details{Used in such problems as Thurstonian scaling.  Although not technically matrix addition, as pointed out by Krus, there are many applications where the sum or difference of two vectors or matrices is a useful operation.  An alternative operation for vectors is  outer(x ,y , FUN="+") but this does not work for matrices.  
}
\value{a n by p matix of sums

}
\references{Krus, D. J. (2001) Matrix addition. Journal of Visual Statistics, 1, (February, 2001).}
\author{William Revelle}

\examples{

x <- seq(1,4)
z <- x \%+\% -t(x)
x
z
#compare with outer(x,-x,FUN="+")
x <- matrix(seq(1,6),ncol=2)
y <- matrix(seq(1,10),nrow=2)
z <- x \%+\% y
x
y
z
#but compare this with outer(x ,y,FUN="+") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
