\name{lowerUpper}
\alias{lowerUpper}

\title{Combine two square matrices to have a lower off diagonal for one, upper off diagonal for the other}

\description{
When reporting correlation matrices for two samples (e.g., males and females), it is convenient to show them as one matrix, with entries below the diagonal representing one matrix, and entries above the diagonal the other matrix.  It is also useful to compare a correlation matrix with the residuals from a fitted (e.g., factor) model.  }

\usage{
lowerUpper(lower, upper=NULL, diff=FALSE)
}

\arguments{
  \item{lower}{A square matrix}
  \item{upper}{A square matrix of the same size as the first (if omitted, then the matrix is converted to two symmetric matrices).}
  \item{diff}{Find the difference between the first and second matrix and put the results in the above the diagonal entries.  }
}
\details{If just one matrix is provided (i.e., upper is missing), it is decomposed into two square matrices, one equal to the lower off diagonal entries, the other to the upper off diagonal entries. In the normal case two symmetric matrices are provided and combined into one non-symmetric matrix with the lower off diagonals representing the lower matrix and the upper off diagonals representing the upper matrix.

If diff is true, the upper off diagonal matrix reflects the differences between the two matrices.
}
\value{Either one matrix or a list of two}

\author{William Revelle}


\seealso{\code{\link{read.clipboard.lower}}, \code{\link{corPlot}}
}
\examples{
 b1 <- Bechtoldt.1
 b2 <- Bechtoldt.2
 b12 <- lowerUpper(b1,b2)
 cor.plot(b12)
 diff12 <- lowerUpper(b1,b2,diff=TRUE)

 corPlot(t(diff12),numbers=TRUE,main="Bechtoldt1 and the differences from Bechtoldt2")
 
 #Compare r and partial r
 lower <- lowerCor(sat.act)
 upper <- partial.r(sat.act)
 both = lowerUpper(lower,upper)
 corPlot(both,numbers=TRUE,main="r and partial r for the sat.act data set")
#now show the differences
 both = lowerUpper(lower,upper,diff=TRUE)
 corPlot(both,numbers=TRUE,main="Differences between r and partial r for the sat.act data set")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }% at least one, from doc/KEYWORDS
