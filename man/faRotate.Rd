\name{faRotations}
\alias{faRotations}
\title{Multiple rotations of factor loadings to find local minima}


\description{A dirty little secret of factor rotation algorithms is the problem of local minima (Nguyen and Waller,2022).  Following ideas in that article, we allow for multiple random restarts and then return the global optimal solution.  Used as part of the \code{\link{fa}} function or available as a stand alone function. 
}
\usage{
faRotations(loadings, r = NULL, rotate = "oblimin", hyper = 0.15, n.rotations = 10,...)
}

\arguments{
  \item{loadings}{Factor loadings matrix from \code{\link{fa}} or \code{\link{pca}} or any N x k loadings matrix}
  \item{r}{The correlation matrix used to find the factors.  (Used to find the factor indeterminancy of the solution)}

 \item{rotate}{"none", "varimax", "quartimax",  "bentlerT", "equamax", "varimin", "geominT" and "bifactor" are orthogonal rotations.  "Promax", "promax", "oblimin", "simplimax", "bentlerQ,  "geominQ" and "biquartimin" and "cluster" are possible oblique transformations of the solution.  Defaults to oblimin. 
 }
  \item{hyper}{The value defining when a loading is in the ``hyperplane".}

  \item{n.rotations}{The number of random restarts to use.}
  \item{...}{additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax}
}


\details{Nguyen and Waller review the problem of local minima in factor analysis.  This is a problem for all rotation algorithms, but is more so for some.  \code{\link{faRotate}} generates n.rotations different starting values and then applies the specified rotation to the original loadings using multiple start values.  Hyperplane counts and complexity indices are reported for each starting matrix, and the one with the highest hyoerplane count and the lowest complexity is returned.
}

\value{
\item{loadings}{The best rotated solution}
\item{Phi}{Factor correlations}
\item{rotation.stats}{Hyperplane count, complexity.}
\item{rot.mat}{The rotation matrix used.}
}
\references{
Nguyen, H. V., & Waller, N. G. (2022, January 6). Local Minima and Factor Rotations in Exploratory Factor Analysis. Psychological Methods. Advance online publication. doi 10.1037/met0000467
}
\author{William Revelle
}
\note{Adapted from the fungible package by Waller
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{fa}} 
}
\examples{
f5 <- fa(psychTools::bfi[,1:25],5,rotate="none")
faRotations(f5,n.rotations=10)   #note that the factor analysis needs to not do the rotation
faRotations(f5$loadings)  #matrix input
geo <- faRotations(f5,rotate="geominQ",n.rotation=10)
 # a popular alternative, but more sensitive to local minima
describe(geo$rotation.stats[,1:3]) 
 }


\keyword{ multivariate }
\keyword{ models}
