\name{omega.graph}
\alias{omega.graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graph hierarchical factor structures }
\description{Hierarchical factor structures represent the correlations between variables in terms of a smaller set of correlated factors which themselves can be represented by a higher order factor.

Two alternative solutions to such structures are found by the omega function.  The correlated factors solutions represents the effect of the higher level, general factor, through its effect on the correlated factors.  The other representation makes use of the Schmid Leiman transformation to find the direct effect of the general factor upon the original variables as well as the effect of orthogonal residual group factors upon the items. 

Graphic presentations of these two alternatives are helpful in understanding the structure.  omega.graph draws both such structures.  Graphs are drawn directly onto the graphics window or expressed in ``dot" commands for conversion to graphics using implementations of Graphviz.
}
\usage{
omega.graph(om.results, out.file = NULL,  sl = TRUE, labels = NULL, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = "RL", digits = 1, title = "Omega", ...)
}

\arguments{
  \item{om.results}{The output from the omega function }
  \item{out.file}{ Optional output file for off line analysis using Graphviz }
  \item{sl}{ Orthogonal clusters using the Schmid-Leiman transform (sl=TRUE) or oblique clusters }
  \item{labels}{ variable labels }
  \item{size}{size of graphics window }
  \item{node.font}{ What font to use for the items}
  \item{edge.font}{What font to use for the edge labels }
  \item{rank.direction}{ Defaults to left to right }
  \item{digits}{ Precision of labels }
  \item{title}{ Figure title }
  \item{\dots}{Other options to pass into the graphics packages }
}
\details{Requires the Rgraphviz package.  omega requires the GPArotation package.
}
\value{
  
  \item{clust.graph }{A graph object}
 }
\references{ \url{http://personality-project.org/r/r.omega.html} \cr

Revelle, W. (1979).  Hierarchical cluster analysis and the internal structure of tests. Multivariate Behavioral Research, 14, 57-74. (\url{http://personality-project.org/revelle/publications/iclust.pdf})

 Zinbarg, R.E., Revelle, W., Yovel, I., & Li. W.  (2005). Cronbach's Alpha, Revelle's Beta, McDonald's Omega: Their relations with each and two alternative conceptualizations of reliability. Psychometrika. 70, 123-133.  \url{http://personality-project.org/revelle/publications/zinbarg.revelle.pmet.05.pdf}

 
 Zinbarg, R., Yovel, I., Revelle, W. & McDonald, R. (2006).  Estimating generalizability to a universe of indicators that all have one attribute in common:  A comparison of estimators for omega.  Applied Psychological Measurement, 30, 121-144. DOI: 10.1177/0146621605278814 \url{http://apm.sagepub.com/cgi/reprint/30/2/121}
}
\author{ 
 \url{http://personality-project.org/revelle.html} \cr
Maintainer: William Revelle  \email{  revelle@northwestern.edu           } 
}
\note{ Requires rgraphviz.   -- omega requires GPARotation
}
\seealso{ \code{\link{omega}}, \code{\link{ICLUST.rgraph}} }
\examples{
om24 <- omega(Harman74.cor$cov,4)  #run omega
om24pn <- omega.graph(om24,sl=FALSE) #show the structure
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(om.results,out.file=NULL,sl=TRUE,labels=NULL,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction="RL", digits=2,title="Omega", ...){
    require(Rgraphviz)
    
   if (sl) {factors <- as.matrix(om.results$schmid$sl)   } else{factors <- as.matrix(om.results$schmid$oblique)}
   rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
 
   num.var <- dim(factors)[1]   #how many variables?
  if (sl) {num.factors <- dim(factors)[2] -3 } else {num.factors <- dim(factors)[2]}
  
   vars <- paste("V",1:num.var,sep="")   
   fact <- c("g",paste("F",1:num.factors,sep=""))
   clust.graph <-  new("graphNEL",nodes=c(vars,fact),edgemode="directed")
   graph.shape <- c(rep("box",num.var),rep("ellipse",num.factors+1))
   graph.rank <- c(rep("sink",num.var),rep("",num.factors+1))
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   edge.label <- rep("",num.var*2)
   edge.name <- rep("",num.var*2)
   names(edge.label) <-  seq(1:num.var*2)  
  #show the cluster structure with ellipses
   if (sl) {
   l <- matrix(factors[,2:(num.factors+1)],ncol=num.factors) } else { l <- factors }
   m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max), 
        ncol = 1)
        
  if (sl) { for (i in 1:num.var) {
     clust.graph <- addEdge(fact[1], vars[i], clust.graph,1) } } else {
        for (i in 1:num.factors) {clust.graph <- addEdge(fact[1], fact[i+1], clust.graph,1) } }
   for (i in 1:num.var) {  clust.graph <- addEdge(fact[1+m1[i]], vars[i], clust.graph,1) }  

 if(FALSE) {      
     edge.label[(i-1)*2+1] <- results[i,"r1"]
     edge.name [(i-1)*2+1]  <- paste(row.names(results)[i],"~", results[i,1],sep="")
     clust.graph <- addEdge(row.names(results)[i], results[i,2], clust.graph,1)
      edge.label[i*2] <- results[i,"r2"]
      edge.name [i*2]  <- paste(row.names(results)[i],"~", results[i,2],sep="")
     } 
     
 nAttrs <- list()  #node attributes
 eAttrs <- list()  #edge attributes


 if (!is.null(labels)) {var.labels <- c(labels,fact)
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
 cluster.vars <- subGraph(fact,clust.graph)
 observed <- list(list(graph=obs.var,cluster=TRUE,attrs=c(rank="")))
 plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed) 
if(!is.null(out.file) ){toDot(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs) }
return(clust.graph)
   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
