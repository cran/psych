\name{superMatrix}
\alias{superMatrix}
\alias{super.matrix}
\title{Form a super matrix from two sub matrices.  }
\description{Given the matrices nXm, and jYk, form the super matrix of dimensions (n+j) and (m+k) with  with elements x and y along the super diagonal. 
Useful when considering structural equations.  The measurement models x and y can be combined into a larger measurement model of all of the variables.  If either x or y is a list of matrices, then recursively form a super matrix of all of those elements.
}
\usage{
superMatrix(x,y)
super.matrix(x, y)  #Deprecated
}

\arguments{
  \item{x}{A n x m matrix or a list of such matrices }
  \item{y}{A j x k matrix or a list of such matrices}
}
\details{Several functions, e.g., \code{\link{sim.structural}},\code{\link{structure.graph}}, \code{\link{make.keys}} use matrices that can be thought of as formed from a set of submatrices.  In particular, when using \code{\link{make.keys}} in order to score a set of items (\code{\link{score.items}}) or to form specified clusters (\code{\link{cluster.cor}}), it is convenient to define different sets of scoring keys for different sets of items and to combine these scoring keys into one super key.

}
\value{
 A (n+j)  x (m +k) matrix with appropriate row and column names
}

\author{William Revelle}

\seealso{ \code{\link{sim.structural}},\code{\link{structure.graph}}, \code{\link{make.keys}}}
\examples{
mx <- matrix(c(.9,.8,.7,rep(0,4),.8,.7,.6),ncol=2)
my <- matrix(c(.6,.5,.4))

colnames(mx) <- paste("X",1:dim(mx)[2],sep="")
rownames(mx) <- paste("Xv",1:dim(mx)[1],sep="")
colnames(my) <- "Y"
rownames(my) <- paste("Yv",1:3,sep="")
mxy <- superMatrix(mx,my)
#show the use of a list to do this as well
key1 <- make.keys(6,list(first=c(1,-2,3),second=4:6,all=1:6))  #make a scoring key
key2 <- make.keys(4,list(EA=c(1,2),TA=c(3,4)))
superMatrix(list(key1,key2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

