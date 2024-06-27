\name{parcels}
\alias{parcels}
\alias{keysort}
\title{Find miniscales (parcels) of size 2 or 3 from a set of items}

\description{Given a set of n items, form n/2 or n/3 mini scales or parcels of the most similar pairs or triplets of items.  These may be used as the basis for subsequent scale construction or multivariate (e.g., factor) analysis.
}
\usage{
parcels(x, size = 3, max = TRUE, flip=TRUE,congruence = FALSE)
keysort(keys)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix/dataframe of items or a correlation/covariance matrix of items}
  \item{size}{Form parcels of size 2 or size 3}
  \item{flip}{if flip=TRUE, negative correlations lead to at least one item being negatively scored}
  \item{max}{Should item correlation/covariance be adjusted for their maximum correlation}
  \item{congruence}{Should the correlations be converted to congruence coefficients?}
  \item{keys}{Sort a matrix of keys to reflect item order as much as possible}
}
\details{Items used in measuring ability or other aspects of personality are typically not very reliable.  One suggestion has been to form items into homogeneous item composites (HICs), Factorially Homogeneous Item Dimensions (FHIDs) or mini scales (parcels).  Parcelling may be done rationally, factorially, or empirically based upon the structure of the correlation/covariance matrix.  \code{link{parcels}} facilitates the finding of parcels by forming a keys matrix suitable for using in \code{\link{score.items}}.  These keys represent the n/2 most similar pairs or the n/3 most similar triplets.

The algorithm is straightforward:  For size = 2, the correlation matrix is searched for the highest correlation.  These two items form the first parcel and are dropped from the matrix.  The procedure is repeated until there are no more pairs to form.

For size=3, the three items with the greatest sum of variances and covariances with each other is found.  This triplet is the first parcel.  All three items are removed and the procedure then identifies the next most similar triplet.  The procedure repeats until n/3 parcels are identified.  
}
\value{
\item{keys}{A matrix of scoring keys to be used to form mini scales (parcels) These will be in order of importance, that is, the first parcel (P1) will reflect the most similar pair or triplet.  The keys may also be sorted by average item order by using the keysort function.}

}
\references{
Cattell, R. B.	(1956).	Validation and intensification of the sixteen personality factor questionnaire.	Journal	of	Clinical	Psychology ,	12 (3),	205	-214.
}
\author{William Revelle}




\seealso{ \code{\link{scoreItems}} to score the parcels or  \code{\link{iclust}} for an alternative way of forming item clusters.
}
\examples{
parcels(Thurstone)
keys <- parcels(bfi)
keys <- keysort(keys)
scoreItems(keys,bfi)
}
\keyword{ multivariate }