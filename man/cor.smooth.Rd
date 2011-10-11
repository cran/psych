\name{cor.smooth}
\alias{cor.smooth}

\title{Smooth a non-positive definite correlation matrix to make it positive definite}

\description{Factor analysis requires positive definite correlation matrices.  Unfortunately, with pairwise deletion of missing data or if using \code{\link{tetrachoric}} or \code{\link{polychoric}} correlations, not all correlation matrices are positive definite.  cor.smooth does a eigenvector (principal components) smoothing.  Negative eigen values are replaced with 100 * .Machine$double.eps, the matrix is reproduced and forced to a correlation matrix using cov2cor.
}
\usage{
cor.smooth(x)
}
\arguments{
  \item{x}{A correlation matrix}
}
\details{The smoothing is done by eigen value decomposition.  eigen values < .Machine$double.eps are changed to 100  * .Machine$double.eps.  The positive eigen values are rescaled to sum to the number of items.  The matrix is recomputed (eigen.vectors \%*\% diag(eigen.values) \%*\% t(eigen.vectors) and forced to a correlation matrix using cov2cor. (See Bock, Gibbons and Muraki, 1988 and Wothke, 1993). 

This does not implement the Knol and ten Berge (1989) solution, nor do nearcor and posdefify in sfmsmisc, not does nearPD in Matrix.  
}
\value{The smoothed matrix with a warning reporting that smoothing was necessary (if smoothing was in fact necessary).
}
\references{
R. Darrell Bock, Robert Gibbons and Eiji Muraki (1988) Full-Information Item Factor Analysis. Applied Psychological Measurement, 12 (3), 261-280.

Werner Wothke (1993), Nonpositive definite matrices in structural modeling. In Kenneth A. Bollen and J. Scott Long (Editors),Testing structural equation models, Sage Publications, Newbury Park. 

D.L. Knol and JMF ten Berge (1989) Least squares approximation of an improper correlation matrix by a proper one.  Psychometrika, 54, 53-61.

}
\author{William Revelle}


\seealso{ \code{\link{tetrachoric}}, \code{\link{polychoric}}, \code{\link{fa}} and \code{\link{irt.fa}}, and the \code{\link{burt}} data set.

See also nearcor and posdefify in the sfsmisc package and \code{\link{nearPD}} in the Matrix package.

}
\examples{
bs <- cor.smooth(burt)  #burt data set is not positive definite
round(burt - bs,3) 
fa(burt,2) #this throws a warning that the matrix yields an improper solution
fa(cor.smooth(burt),2)  #This throws a warning that the matrix was improper, but produces a better solution 

}

\keyword{ multivariate }
\keyword{ models}
