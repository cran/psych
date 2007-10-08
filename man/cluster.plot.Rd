\name{cluster.plot}
\alias{cluster.plot}

\title{Plot factor/cluster loadings and assign items to clusters by their highest loading}
\description{Cluster analysis and factor analysis are procedures for grouping items in terms of a smaller number of (latent) factors or (observed) clusters.  Graphical presentations of clusters typically show tree structures, although they can be represented in terms of item by cluster correlations.  

Cluster.plot plots items by their cluster loadings (taken, e.g., from ICLUST).  Cluster membership may be assigned apriori or may be determined in terms of the highest (absolute) cluster loading for each item.  
}
\usage{
cluster.plot(ic.results, cluster = NULL, cut = 0, labels=NULL,title = "Cluster plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ic.results}{A factor analysis or cluster analysis output including the loadings, or a matrix of item by cluster correlations }
  \item{cluster}{ A vector of cluster membership }
 
  \item{cut}{ Assign items to clusters if the absolute loadings are > cut }
  \item{labels}{If row.names exist they will be added to the plot, or, if they don't, labels can be specified.  If labels =NULL, and there are no row names, then variables are labeled by row number.)}
  \item{title}{ Any title}
  
}
\details{
  
}
\value{
 Graphical output is presented
}

\author{William Revelle}
\note{ September 5, 2007

 
}
\seealso{ \code{\link{ICLUST}},\code{\link{ICLUST.graph}}}
\examples{circ.data <- circ.sim(24,500)
circ.fa <- factor.pa(circ.data,2)
cluster.plot(circ.fa,cut=.5)
}

\keyword{multivariate }
\keyword{hplot }% __ONLY ONE__ keyword per line
\keyword{cluster}
