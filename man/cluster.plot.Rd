\name{cluster.plot}
\alias{cluster.plot}
\alias{fa.plot}
\alias{factor.plot}

\title{Plot factor/cluster loadings and assign items to clusters by their highest loading.}
\description{Cluster analysis and factor analysis are procedures for grouping items in terms of a smaller number of (latent) factors or (observed) clusters.  Graphical presentations of clusters typically show tree structures, although they can be represented in terms of item by cluster correlations.  

Cluster.plot plots items by their cluster loadings (taken, e.g., from \code{\link{ICLUST}}) or factor loadings (taken, eg., from \code{\link{fa}}).  Cluster membership may be assigned apriori or may be determined in terms of the highest (absolute) cluster loading for each item.  

If the input is an object of class "kmeans", then the cluster centers are plotted. 
}
\usage{


cluster.plot(ic.results, cluster = NULL, cut = 0, labels=NULL,
          title = "Cluster plot",pch=18,pos,...)
fa.plot(ic.results, cluster = NULL, cut = 0, labels=NULL,title, 
             jiggle=FALSE,amount=.02,pch=18,pos,...)
factor.plot(ic.results, cluster = NULL, cut = 0, labels=NULL,title,jiggle=FALSE,
                  amount=.02,pch=18,pos,...)  #deprecated
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ic.results}{A factor analysis or cluster analysis output including the loadings, or a matrix of item by cluster correlations. Or the output from a kmeans cluster analysis. }
  \item{cluster}{ A vector of cluster membership }
  
  \item{cut}{ Assign items to clusters if the absolute loadings are > cut }
  \item{labels}{If row.names exist they will be added to the plot, or, if they don't, labels can be specified.  If labels =NULL, and there are no row names, then variables are labeled by row number.)}
  \item{title}{ Any title}
  \item{jiggle}{When plotting with factor loadings that are almost identical, it is sometimes useful to "jiggle" the points by jittering them. The default is to not jiggle.}
  \item{amount}{if jiggle=TRUE, then how much should the points be jittered?}
  \item{pch}{factor and clusters are shown with different pch values, starting at pch+1}
  \item{pos}{Position of the text for labels for two dimensional plots. 1=below, 2 = left, 3 = above, 4= right}
  \item{...}{Further options to plot}
  
}
\details{Results of either a factor analysis or cluster analysis are plotted.  Each item is assigned to its highest loading factor, and then identified by variable name as well as cluster (by color). 
The cluster assignments can be specified to override the automatic clustering by loading.
Both of these functions may be called directly or by calling the generic plot function.  (see example).
  
}
\value{
 Graphical output is presented.
}

\author{William Revelle}

\seealso{ \code{\link{ICLUST}}, \code{\link{ICLUST.graph}}, \code{\link{fa.graph}}, \code{\link{plot.psych}}}
\examples{
circ.data <- circ.sim(24,500)
circ.fa <- fa(circ.data,2)
plot(circ.fa,cut=.5)
}

\keyword{multivariate }
\keyword{hplot }% __ONLY ONE__ keyword per line
\keyword{cluster}
