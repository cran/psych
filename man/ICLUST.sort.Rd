\name{ICLUST.sort}
\alias{ICLUST.sort}
\title{Sort items by absolute size of cluster loadings}
\description{Given a cluster analysis or factor analysis loadings matrix, sort the items by the (absolute) size of each column of loadings.  Used as part of ICLUST and SAPA analyses.
}
\usage{
ICLUST.sort(ic.load, cut = 0, labels = NULL,keys=FALSE)
}

\arguments{
  \item{ic.load}{ The output from a factor or principal components analysis, or from ICLUST, or a matrix of loadings.}
  \item{cut}{Do not include items in clusters with absolute loadings less than cut}
  \item{labels}{labels for each item.}
 
  \item{keys}{should cluster keys be returned? Useful if clusters scales are to be scored.}
}
\details{When interpreting cluster or factor analysis outputs, is is useful to group the items in terms of which items have their biggest loading on each factor/cluster and then to sort the items by size of the absolute factor loading.

A stable cluster solution will be one in which the output of these cluster definitions does not vary when clusters are formed from the clusters so defined.

With the keys=TRUE option, the resulting cluster keys may be used to score the original data or the correlation matrix to form clusters from the factors.
}
\value{
  \item{sorted }{A data.frame of item numbers, item contents, and item x factor loadings.}
  \item{cluster}{A matrix of -1, 0, 1s defining each item by the factor/cluster with the row wise largest absolute loading. }

  ...
}
\references{ \url{http://personality-project.org/r/r.ICLUST.html} }
\author{William Revelle }

\note{ Although part of the ICLUST set of programs, this is also more useful for factor or principal components analysis.

}
\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{factor2cluster} }}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
