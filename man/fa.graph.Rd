\name{fa.graph}
\alias{fa.graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graph factor loading matrices}
\description{Factor analysis or principal components analysis results are typically interpreted in terms of the major loadings on each factor.  These structures may be represented as a table of loadings or graphically, where all loadings with an absolute value > some cut point are represented as an edge (path).
}
\usage{
fa.graph(fa.results, out.file = NULL, labels = NULL, cut = 0.3, simple = TRUE, size = c(8, 6), node.font = c("Helvetica", 14), edge.font = c("Helvetica", 10), rank.direction = "RL", digits = 1, title = "Factor Analysis", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fa.results}{The output of factor analysis or principal components analysis }
  \item{out.file}{ If it exists, a dot representation of the graph will be stored here }
  \item{labels}{ Variable labels }
  \item{cut}{ Loadings with abs(loading) > cut will be shown }
  \item{simple}{Only the biggest loading per item is shown}
  \item{size}{  }
  \item{node.font}{  }
  \item{edge.font}{ ~ }
  \item{rank.direction}{  }
  \item{digits}{ Number of digits to show as an edgelable }
  \item{title}{ Graphic title }
  \item{\dots}{ other parameters }
}
\details{Path diagram representations have become standard in confirmatory factor analysis, but are not yet common in exploratory factor analysis.  Representing factor structures graphically helps some people understand the structure. 
}
\value{ A graph is drawn using rgraphviz.  If an output file is specified, the graph instructions are also saved in the dot language. }
  


\author{William Revelle }
\note{ Requires rgraphviz}

 
}
\seealso{  \code{\link{omega.graph}}, \code{\link{ICLUST.graph}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
# test.simple <- factor.pa(item.sim(16),2)
#fa.graph(test.simple)

## The function is currently defined as
function(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction="RL", digits=1,title="Factor Analysis", ...){
    require(Rgraphviz)
    
factors <- as.matrix(fa.results$loadings)
   rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
  
   num.var <- dim(factors)[1]   #how many variables?
   if (is.null(num.var) ){num.var <- length(factors)
       num.factors <- 1} else {
   num.factors <- dim(factors)[2]}
   if (simple) {k=1} else {k <- num.factors}
   vars <- paste("V",1:num.var,sep="")   
   fact <- paste("F",1:num.factors,sep="")
   clust.graph <-  new("graphNEL",nodes=c(vars,fact),edgemode="directed")
   graph.shape <- c(rep("box",num.var),rep("ellipse",num.factors))
   graph.rank <- c(rep("sink",num.var),rep("",num.factors))
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   edge.label <- rep("",num.var*k)
   edge.name <- rep("",num.var*k)
   names(edge.label) <-  seq(1:num.var*k) 
            
  #show the cluster structure with ellipses

  l <- factors 
  if (num.factors ==1) { 
       
    for (i in 1:num.var) { clust.graph <- addEdge(fact[1], vars[i], clust.graph,1) 
                         edge.label[i] <- round(factors[i],digits)
                         edge.name[i] <- paste(fact[1],"~",vars[i],sep="")
                        }  
      } else { 
     if(simple){   #very simple structure is one loading per variable
        m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max), 
        ncol = 1)       
        for (i in 1:num.var) {clust.graph <- addEdge(fact[m1[i]], vars[i], clust.graph,1) 
                         edge.label[i] <- round(factors[i,m1[i]],digits)
                         edge.name[i] <- paste(fact[m1[i]],"~",vars[i],sep="")
                        }  
                  }   else {   #all loadings > cut in absolute value
                   k <- 1
                   for (i in 1:num.var) {
                   for (f in 1:num.factors) { if (abs(factors[i,f]) > cut) {clust.graph <- addEdge(fact[f], vars[i], clust.graph,1) 
                         edge.label[k] <- round(factors[i,f],digits)
                         edge.name[k] <- paste(fact[f],"~",vars[i],sep="")
                         k <- k+1 }   #end of if 
                         }  #end of factor
                        } # end of variable loop
            }  #end of if simple else
       }   #end of if num.factors ==1 
                  
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
 attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE),graph=list(rankdir=rank.direction, fontsize=10,bgcolor="white" ))
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

