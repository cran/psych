\name{ICLUST.rgraph}
\alias{ICLUST.rgraph}

\title{ Draw an ICLUST graph using the Rgraphviz package }
\description{Given a cluster structure determined by \code{\link{ICLUST}}, create a rgraphic directly using Rgraphviz.  To create dot code to describe the \code{\link{ICLUST}} output with more precision, use \code{\link{ICLUST.graph}}.  As an option, dot code is also generated and saved in a file. To  use the dot code, use either \url{http://www.graphviz.org/} Graphviz or a commercial viewer (e.g., OmniGraffle).
}
\usage{
ICLUST.rgraph(ic.results, out.file = NULL, min.size = 1, short = FALSE, labels = NULL, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = "RL", digits = 2, title = "ICLUST", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
   \item{ic.results}{output list from ICLUST }
   \item{out.file}{ File name to save optional dot code. }
  \item{min.size}{draw a smaller node (without all the information) for clusters < min.size -- useful for large problems}
  \item{short}{if short==TRUE, don't use variable names}
  \item{labels}{vector of text labels (contents) for the variables}
  \item{size}{size of output }
  \item{node.font}{ Font to use for nodes in the graph }
  \item{edge.font}{ Font to use for the labels of the arrows (edges)}
  \item{rank.direction}{LR or RL  }
  \item{digits}{ number of digits to show }
  \item{title}{ any title }
  \item{\dots}{ other options to pass }
}
\details{
  Will create (or overwrite) an output file and print out the dot code to show a cluster structure. This dot file may be imported directly into a dot viewer (e.g.,  \url{http://www.graphviz.org/}).  The "dot" language is a powerful graphic description language that is particulary appropriate for viewing cluster output.  Commercial graphics programs (e.g., OmniGraffle) can also read (and clean up) dot files.  
  
ICLUST.graph takes the output from \code{\link{ICLUST}} results and processes it to provide a pretty picture of the results.  Original variables shown as rectangles and ordered on the left hand side (if rank direction is RL) of the graph.  Clusters are drawn as ellipses and include the alpha, beta, and size of the cluster.  Edges show the cluster intercorrelations.

It is possible to trim the output to not show all cluster information. Clusters < min.size are shown as small ovals without alpha, beta, and size information.
}
\value{Output is a set of dot commands written either to console or to the output file.  These commands may then be used as input to any "dot" viewer, e.g., Graphviz.

Additional output is drawn to main graphics screen. }
\references{ ICLUST: \url{http://personality-project.org/r/r.iclust.html}}
\author{
 \email{revelle@northwestern.edu } \cr
 \url{http://personality-project.org/revelle.html}}

\seealso{  \code{\link{VSS.plot}},  \code{\link{ICLUST}}}
\note{ Requires Rgraphviz}
}

\examples{
test.data <- Harman74.cor$cov
ic.out <- ICLUST(test.data)
ICLUST.rgraph(ic.out)   

## The function is currently defined as
function(ic.results,out.file = NULL, min.size=1,short=FALSE,labels=NULL,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction="RL", digits=2,title="ICLUST", ...){
    require(Rgraphviz)
   clusters <- as.matrix(ic.results$clusters)  
   results <- ic.results$results 
  
   rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
 
 
   #create the items as boxes  
   #add the sign from the clusters 
   num.var <- dim(clusters)[1]   #how many variables?
   num.clust <- num.var - dim(clusters)[2] 
  
   vars <- paste("V",1:num.var,sep="")   
   clust <- paste("C",1:num.clust,sep="")
   clust.graph <-  new("graphNEL",nodes=c(vars,clust),edgemode="directed")
   graph.shape <- c(rep("box",num.var),rep("ellipse",num.clust))
   graph.rank <- c(rep("sink",num.var),rep("",num.clust))
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   edge.label <- rep("",num.clust*2)
   edge.name <- rep("",num.clust*2)
   names(edge.label) <-  seq(1:num.clust*2)  
  #show the cluster structure with ellipses
  for (i in 1:num.clust) {if(results[i,1]>0) { #avoid printing null results
     clust.graph <- addEdge(row.names(results)[i], results[i,1], clust.graph,1)
     edge.label[(i-1)*2+1] <- results[i,"r1"]
     edge.name [(i-1)*2+1]  <- paste(row.names(results)[i],"~", results[i,1],sep="")
     clust.graph <- addEdge(row.names(results)[i], results[i,2], clust.graph,1)
      edge.label[i*2] <- results[i,"r2"]
      edge.name [i*2]  <- paste(row.names(results)[i],"~", results[i,2],sep="")
     }}
 nAttrs <- list()  #node attributes
 eAttrs <- list()  #edge attributes

 if (!is.null(labels)) {var.labels <- c(labels,row.names(results)) #note how this combines variable labels with the cluster variables
  names(var.labels) <-  nodes(clust.graph)
  nAttrs$label <- var.labels
  names(edge.label) <- edge.name
  } 
    names(edge.label) <- edge.name
 nAttrs$shape <- graph.shape
 nAttrs$rank <- graph.rank
 eAttrs$label <- edge.label
 attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE),graph=list(rankdir="RL", fontsize=10,bgcolor="white" ))
 obs.var <- subGraph(vars,clust.graph)
 cluster.vars <- subGraph(clust,clust.graph)
 observed <- list(list(graph=obs.var,cluster=TRUE,attrs=c(rank="sink")))
 plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed) 
 if (!is.null(out.file)) {toDot(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed) }
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate}
\keyword{ cluster}% __ONLY ONE__ keyword per line
