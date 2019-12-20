\name{Pinv}
\alias{Pinv}
\title{Compute the  Moore-Penrose Pseudo Inverse of a matrix}

\description{Given a matrix of less than full rank, the conventional inverse function will fail.  The pseudoinverse or generalized inverse resolves this problem by using just the postive  values of the singular value decomposition d matrix. An adaptation of the ginv function from MASS and the pinv function from pracma. 
}
\usage{
Pinv(X, tol = sqrt(.Machine$double.eps))
}

\arguments{
\item{X}{A correlation or covariance matrix to analyze}
\item{tol}{A very small number. Reject values with eigen values less than tolerance}
}

\details{ The singular value decomposition of a matrix X is UdV where for full rank matrices, d is the vector of eigen values and U and V are the matrices of eigen vectors. The inverse is just U/d.  If the matrix is less than full rank, many of the d values are effectively zero (at the limit of computational accuracy.) Thus, to solve matrix equations with matrices of less than full rank (e.g. the \code{\link{schmid}} Schmid-Leiman solution), we need to find the generalized inverse. 



}
\value{The generalized inverse}

\references{
Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition. Springer. p.100.
}
\author{
William Revelle
}
\note{Adapted from the ginv function in MASS and the pinv function in pracma. Installed here to avoid loading those packages.
}


\seealso{
\code{\link{schmid}}, \code{\link{faCor}}
}
\examples{
round(Pinv(Thurstone) \%*\% Thurstone,2)  #an identity matrix
sl <- schmid(Thurstone,3)  #The schmid-leiman solution is less than full rank
F <- sl$sl[,1:4]    #the SL solution is general + 3 gropus
R <- Thurstone      #
diag(R) <- sl$sl[,5]  #the reproduced matrix (R - U2)
S <- t(Pinv(t(F) \%*\% F) \%*\% t(F) \%*\% R)  #the structure matrix
Phi <- t(S) \%*\%  F \%*\% Pinv(t(F) \%*\% F)  #the factor covariances

}

\keyword{multivariate}