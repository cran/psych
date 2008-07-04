\name{cluster2keys}
\alias{cluster2keys}

\title{Convert a cluster vector (from e.g., kmeans) to a keys matrix suitable for scoring item clusters. }
\description{The output of the kmeans clustering function produces a vector of cluster membership.  The \code{\link{score.items}} and \code{\link{cluster.cor}} functions require a matrix of keys.  cluster2keys does this.
}
\usage{
cluster2keys(c)
}
\arguments{
  \item{c}{A vector of cluster assignments or an object of class ``kmeans" that contains a vector of clusters. }
}
\details{Note that because kmeans will not reverse score items, the clusters defined by kmeans will not necessarily match those of ICLUST with the same number of clusters extracted. 
}
\value{
 
  \item{keys}{A matrix of keys suitable for score.items or cluster.cor}
 
}

\author{William Revelle}

\seealso{ \code{\link{cluster.cor}},\code{\link{score.items}}}
\examples{
test.data <- Harman74.cor$cov
kc <- kmeans(test.data,4)
keys <- cluster2keys(kc)
keys  #these match those found by ICLUST
cluster.cor(keys,test.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }

