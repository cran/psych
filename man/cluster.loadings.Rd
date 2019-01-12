\name{cluster.loadings}
\alias{cluster.loadings}
\title{ Find item by cluster correlations, corrected for overlap and reliability }
\description{
 Given a n x n correlation matrix and a n x c matrix of -1,0,1 cluster weights for those n items on  c clusters, find the correlation of each item with each cluster.  If the item is part of the cluster, correct for item overlap.  Part of the \code{\link{ICLUST}} set of functions, but useful for many item analysis problems.
}
\usage{
cluster.loadings(keys, r.mat, correct = TRUE,SMC=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keys}{ Cluster keys: a matrix of -1,0,1 cluster weights}
  \item{r.mat}{ A correlation matrix }
  \item{correct}{Correct for reliability}
  \item{SMC}{Use the squared multiple correlation as a communality estimate, otherwise use the greatest correlation for each variable}
}
\details{Given a set of items to be scored as (perhaps overlapping) clusters and the intercorrelation matrix of the items, find the clusters and then the correlations of each item with each cluster.  Correct for item overlap by replacing the item variance with its average within cluster inter-item correlation.  

Although part of ICLUST, this may be used in any SAPA (\url{https://sapa-project.org}) application where we are interested in item-whole correlations of items and composite scales.

For information about SAPA see Revelle et al, 2010, 2016.  For information about SAPA based measures of ability, see \url{https://icar-project.com}.

These loadings are particularly interpretable when sorted by absolute magnitude for each cluster (see \code{\link{ICLUST.sort}}). 
}
\value{
  \item{loadings }{A matrix of item-cluster correlations (loadings)}
  \item{cor }{Correlation matrix of the clusters}
  \item{corrected }{Correlation matrix of the clusters, raw correlations below the diagonal, alpha on diagonal, corrected for reliability above the diagonal}
  \item{sd }{Cluster standard deviations}
  \item{alpha }{alpha reliabilities of the clusters}
  \item{G6}{G6* Modified estimated of Guttman Lambda 6}
  \item{count}{Number of items in the cluster}
}
\references{ ICLUST: \url{https://personality-project.org/r/r.ICLUST.html}

Revelle, W., Wilt, J.,  and Rosenthal, A. (2010)  Individual Differences in Cognition: New Methods for examining the Personality-Cognition Link In Gruszka, A.  and Matthews, G. and Szymura, B. (Eds.) Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control, Springer.

 Revelle, W,  Condon, D.M.,  Wilt, J.,  French, J.A., Brown, A.,  and  Elleman, L.G. (2016) Web and phone based data collection using planned missing designs. In  Fielding, N.G.,  Lee, R.M. and  Blank, G. (Eds). SAGE Handbook of Online Research Methods (2nd Ed), Sage Publcations.
 }
\author{Maintainer: William Revelle \email{revelle@northwestern.edu}
   }
        
\note{ Although part of ICLUST, this may be used in any SAPA application where we are interested in item- whole correlations of items and composite scales.} 

\seealso{ \code{\link{ICLUST}}, \code{\link{factor2cluster}},  \code{\link{cluster.cor}} }
\examples{

 r.mat<- Harman74.cor$cov
 clusters <- matrix(c(1,1,1,rep(0,24),1,1,1,1,rep(0,17)),ncol=2)
 cluster.loadings(clusters,r.mat)

 

}
\keyword{multivariate }% at least one, from doc/KEYWORDS
\keyword{ cluster }% __ONLY ONE__ keyword per line
