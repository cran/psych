\name{super.matrix}
\alias{super.matrix}
\title{Form a super matrix from two sub matrices.  }
\description{Given the matrices nXm, and jYk, form the super matrix of dimensions (n+j) and (m+k) with  with elements x and y along the super diagonal. 
Useful when considering structural equations.  The measurement models x and y can be combined into a larger measurement model of all of the variables.
}
\usage{
super.matrix(x, y)
}

\arguments{
  \item{x}{A n x m matrix }
  \item{y}{A j x k matrix }
}

\value{
 A (n+j)  x (m +k) matrix with appropriate row and column names
}

\author{William Revelle}

\seealso{ \code{\link{sim.structural}},\code{\link{structure.graph}}}
\examples{
mx <- matrix(c(.9,.8,.7,rep(0,4),.8,.7,.6),ncol=2)
my <- matrix(c(.6,.5,.4))
colnames(mx) <- paste("X",1:dim(mx)[2],sep="")
rownames(mx) <- paste("Xv",1:dim(mx)[1],sep="")
colnames(my) <- "Y"
rownames(my) <- paste("Yv",1:3,sep="")
super.matrix(mx,my)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

