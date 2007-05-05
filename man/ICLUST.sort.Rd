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
\author{William Revelle \cr
        Department of Psychology \cr
        Northwestern University \cr
        Evanston, Illinois \cr
        \email{  revelle@northwestern.edu           }\cr
 \url{http://personality-project.org/revelle.html}   }
\note{ Although part of the ICLUST set of programs, this is generally more useful for factor or principal components analysis.

}
\seealso{ \code{\link{ICLUST.graph}},\code{\link{ICLUST.cluster}}, \code{\link{cluster.fit} }, \code{\link{VSS}}, \code{\link{factor2cluster} }}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ic.load,cut=0,labels=NULL,loading=FALSE) {
    if (length(labels)==0) {
    var.labels <- rownames(ic.load)} else {var.labels=labels}
    if(loading) {loadings <- ic.load$loadings} else {loadings <- ic.load}
     nclust <- dim(loadings)[2]
     nitems <- dim(loadings)[1]
   
  load <- data.frame(item=seq(1:nitems),content=var.labels,cluster=rep(0,nitems),loadings)
  
  #first find the maximum for each row and assign it to that cluster
   load$cluster <- apply(abs(loadings),1,which.max)
  for (i in 1:nitems) {if (abs(loadings[i,load$cluster[i]]) < cut) {load$cluster[i] <-nclust+1}}
 
  ord <- sort(load$cluster,index.return=TRUE)
  load[1:nitems,] <- load[ord$ix,]
  rownames(load)[1:nitems] <- rownames(load)[ord$ix]
  
  items <- c(table(load$cluster),1)   #how many items are in each cluster?
  #now sort the loadings that have their highest loading on each cluster
   first <- 1
	for (i in 1:nclust) {
	last <- first + items[i]- 1
	ord <- sort(abs(load[first:last,i+3]),decreasing=TRUE,index.return=TRUE)
   load[first:last,] <- load[ord$ix+first-1,]
    rownames(load)[first:last] <- rownames(load)[ord$ix+first-1]
    first <- first + items[i]
    }
    if (first < nitems) load[first:nitems,"cluster"] <- 0   #assign items less than cut to 0
 ICLUST.sort <- load }
  
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
