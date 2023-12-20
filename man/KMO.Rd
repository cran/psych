\name{KMO}
\alias{KMO}
\title{Find the Kaiser, Meyer, Olkin Measure of Sampling Adequacy}

\description{Henry Kaiser (1970) introduced an Measure of Sampling Adequacy (MSA) of factor analytic data matrices. Kaiser and Rice (1974) then modified it. This is just a function of the squared elements of the `image' matrix compared to the squares of the original correlations.  The overall MSA as well as estimates for each item are found. The index is known as the Kaiser-Meyer-Olkin (KMO) index.}
\usage{
KMO(r)
}
\arguments{
  \item{r}{A correlation matrix or a data matrix (correlations will be found)}
}

\details{Let  \eqn{S^2 = diag(R^{-1})^{-1} } and \eqn{Q = SR^{-1}S}.  Then Q is said to  be the anti-image intercorrelation matrix.  Let \eqn{sumr2 = \sum{R^2}} and \eqn{sumq2 = \sum{Q^2}} for all off diagonal elements of R and Q, then  \eqn{SMA=sumr2/(sumr2 + sumq2)}.  Although originally MSA was 1 - sumq2/sumr2  (Kaiser, 1970), this was modified in Kaiser and Rice, (1974) to be   \eqn{SMA=sumr2/(sumr2 + sumq2)}.  This is the formula used by Dziuban and Shirkey (1974) and by SPSS.
	
In his delightfully flamboyant style, Kaiser (1975)
suggested that KMO > .9 were marvelous, in the .80s, mertitourious, in the .70s, middling, in the .60s, medicore, in the 50s, miserable, and  less than .5, unacceptable.
	
	An alternative measure of whether the matrix is factorable is the Bartlett test  \code{\link{cortest.bartlett}} which tests the degree that the matrix deviates from an identity matrix.
	
	}
\value{
\itemize{
\item MSA: The overall Measure of Sampling Adequacy
\item MSAi:  The measure of sampling adequacy for each item
\item Image: The Image correlation matrix (Q)
}
}
\references{
H.~F. Kaiser. (1970) A second generation little jiffy.
Psychometrika, 35(4):401--415.

H.~F. Kaiser and J.~Rice. (1974) Little jiffy, mark iv.
Educational and Psychological Measurement, 34(1):111--117.

H.F. Kaiser. 1974) An index of factor simplicity.  Psychometrika,  39 (1) 31-36.

Dziuban, Charles D. and Shirkey, Edwin C. (1974) When is a correlation matrix appropriate for factor analysis? Some decision rules. Psychological Bulletin, 81 (6) 358 - 361.

}
\author{
William Revelle
}


\seealso{ See Also as \code{\link{fa}},  \code{\link{cortest.bartlett}}, \code{\link{Harman.political}}.
}
\examples{
KMO(Thurstone)
KMO(Harman.political)   #compare to the results in Dziuban and Shirkey (1974)

}

\keyword{ multivariate }
\keyword{ models}
