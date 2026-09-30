\name{Pinv}
\alias{Pinv}
\alias{matSqrt}
\alias{invMatSqrt}
\title{Compute the  Moore-Penrose Pseudo Inverse of a matrix and other tools}

\description{Given a matrix of less than full rank, the conventional inverse function will fail.  The pseudoinverse or generalized inverse resolves this problem by using just the postive  values of the singular value decomposition d matrix. This is an adaptation of the ginv function from MASS and the pinv function from pracma.   Includes several other simple matrix tools such as matSqrt and invMatSqrt which can be used for ``whitening" a matrix.
}
\usage{
Pinv(R, tol = sqrt(.Machine$double.eps))
matSqrt(x) #take the sqrt of a matrix
invMatSqrt(x) #find the sqrt of the inverse matrix

}

\arguments{
\item{R}{A correlation or covariance matrix to analyze}
\item{x}{A correlation or covariance matrix to analyze}
\item{tol}{A very small number. Reject values with eigen values less than tolerance}
}

\details{ The singular value decomposition of a matrix X is UdV where for full rank matrices, d is the vector of eigen values and U and V are the matrices of eigen vectors. The inverse is just U(1/d) U.  If the matrix is less than full rank, many of the d values are effectively zero (at the limit of computational accuracy.) Thus, to solve matrix equations with matrices of less than full rank (e.g. the \code{\link{schmid}} Schmid-Leiman solution), we need to find the generalized inverse. 

\code{\link{matSqrt}} and \code{\link{invMatSqrt}} will find the square root of a matrix or the square root of the \code{\link{Pinv}} (Penrose pseudo inverse) of a matrix.

\code{\link{invMatSqrt}} can be used to find the ``whitening" matrix of a correlation or covariances such that for a correlation matrix, R, W= invMatSqrt(R): 
\deqn{W = R^{-.5}}  
and 
\deqn{W R W' = I.}  
    

}
\value{The generalized inverse,  or the sqrt of the generalized inverse, or the sqrt of the  Penrose inverse.}

\references{

Agnan Kessy, Alex Lewin and  Korbinian Strimmer (2018) Optimal Whitening and Decorrelation, The American Statistician, 72:4, 309-314.

Marco Del Giudice (2024) Data matrix disattenuation: a simple, effective  method for correcting measurement error in multivariate datasets.  Preprint: PsyArXiv 


Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition. Springer. p.100.Optimal

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


W <- invMatSqrt(Thurstone) #find the whitening of Thurstone
round(W \%*\% Thurstone \%*\% W,2)

}

\keyword{multivariate}