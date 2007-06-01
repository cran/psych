\name{ICLUST.sort}
\alias{ICLUST.sort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sort items by absolute size of cluster loadins}
\description{Given a cluster analysis or factor analysis loadings matrix, sort the items by the (absolute) size of each column of loadings.  Used as part of ICLUST and SAPA analyses.
}
\usage{
ICLUST.sort(ic.load, cut = 0, labels = NULL,loading=TRUE)
}

\arguments{
  \item{ic.load}{ A loading matrix from a factor or principal components analysis, or from ICLUST. }
  \item{cut}{Do not include items in clusters with absolute loadings less than cut}
  \item{labels}{labels for each item.}
  \item{loading}{if coming from a factor analysis or ICLUST output, use loading=TRUE, but, if just sorting a matrix use loading=FALSE}
}
\details{When interpreting cluster or factor analysis outputs, is is useful to group the items in terms of which items have their biggest loading on each factor/cluster and then to sort the items by size of the absolute factor loading.

A stable cluster solution will be one in which the output of these cluster definitions does not vary when clusters are formed from the clusters so defined.
}
\value{
  \item{cluster }{A matrix of -1, 0, 1s defining each item by the factor/cluster with the row wise largest absolute loading. }
  \item{load }{A data.frame of item numbers, item contents, and item x factor loadings.}
  ...
}
\references{ \url{http://personality-project.org/r/r.ICLUST.html} }
\author{William Revelle }

\note{ Although part of the ICLUST set of programs, this is generally more useful for factor or principal components analysis.

}
\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{factor2cluster} }}
\examples{

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
