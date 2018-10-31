\name{cor.smooth}
\alias{cor.smooth}
\alias{cor.smoother}

\title{Smooth a non-positive definite correlation matrix to make it positive definite}

\description{Factor analysis requires positive definite correlation matrices.  Unfortunately, with pairwise deletion of missing data or if using \code{\link{tetrachoric}} or \code{\link{polychoric}} correlations, not all correlation matrices are positive definite.  cor.smooth does a eigenvector (principal components) smoothing.  Negative eigen values are replaced with 100  * eig.tol, the matrix is reproduced and forced to a correlation matrix using cov2cor.
}
\usage{
cor.smooth(x,eig.tol=10^-12)
cor.smoother(x,cut=.01)
}
\arguments{
  \item{x}{A correlation matrix or a raw data matrix.}
  \item{eig.tol}{the minimum acceptable eigenvalue}.
  \item{cut}{Report all abs(residuals)  > cut}
}
\details{The smoothing is done by eigen value decomposition.  eigen values < eig.tol are changed to 100  * eig.tol.  The positive eigen values are rescaled to sum to the number of items.  The matrix is recomputed (eigen.vectors \%*\% diag(eigen.values) \%*\% t(eigen.vectors) and forced to a correlation matrix using cov2cor. (See Bock, Gibbons and Muraki, 1988 and Wothke, 1993). 

This does not implement the Knol and ten Berge (1989) solution, nor do nearcor and posdefify in sfmsmisc, not does nearPD in Matrix. As Martin Maechler puts it in the posdedify function, "there are more sophisticated algorithms to solve this and related problems."  

cor.smoother examines all of nvar minors of rank nvar-1 by systematically dropping one variable at a time and finding the eigen value decomposition.  It reports those variables, which, when dropped, produce a positive definite matrix.  It also reports the number of negative eigenvalues when each variable is dropped.  Finally, it compares the original correlation matrix to the smoothed correlation matrix and reports those items with absolute deviations great than cut.  These are all hints as to what might be wrong with a correlation matrix.  
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

See also nearcor and posdefify in the sfsmisc package and nearPD in the Matrix package.

}
\examples{
bs <- cor.smooth(burt)  #burt data set is not positive definite
plot(burt[lower.tri(burt)],bs[lower.tri(bs)],ylab="smoothed values",xlab="original values")
abline(0,1,lty="dashed")

round(burt - bs,3) 
fa(burt,2) #this throws a warning that the matrix yields an improper solution
#Smoothing first throws a warning that the matrix was improper, 
#but produces a better solution 
fa(cor.smooth(burt),2)  

#this next example is a correlation matrix from DeLeuw used as an example 
#in Knol and ten Berge.  
#the example is also used in the nearcor documentation
 cat("pr is the example matrix used in Knol DL, ten Berge (1989)\n")
 pr <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
		0.477, 1,     0.516, 0.233, 0.682, 0.75,
		0.644, 0.516, 1,     0.599, 0.581, 0.742,
		0.478, 0.233, 0.599, 1,     0.741, 0.8,
		0.651, 0.682, 0.581, 0.741, 1,     0.798,
		0.826, 0.75,  0.742, 0.8,   0.798, 1),
	      nrow = 6, ncol = 6)
	      
sm <- cor.smooth(pr)
resid <- pr - sm
# several goodness of fit tests
# from Knol and ten Berge
tr(resid \%*\% t(resid)) /2

# from nearPD
sum(resid^2)/2

}

\keyword{ multivariate }
\keyword{ models}
