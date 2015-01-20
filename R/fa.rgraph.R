#Created May 24, 2007
#modifed April 12, 2008 to allow figures from matrices that are not loadings
#take the output from a factor analysis and graph it using rgraphviz
#Except for the ability to write a dot file, this has been replaced by fa.diagram to avoid using Rgraphviz (September, 2009)
"fa.rgraph" <-
function(fa.results,out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction=c("RL","TB","LR","BT"), digits=1,main="Factor Analysis",graphviz=TRUE, ...){
    	if (!requireNamespace('Rgraphviz')) {stop("I am sorry, you need to have loaded the Rgraphviz package")
    	 #create several dummy functions to get around the "no visible global function definition" problem
    	nodes <- function() {}
    	addEdge <- function() {}
    	subGraph <- function(){} }
    
  Phi <- NULL  #the default case
 if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {factors <- as.matrix(fa.results$loadings)
                if(!is.null(fa.results$Phi)) Phi <- fa.results$Phi} else {factors <- fa.results}
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
   graph.rank <- c(rep("sink",num.var),rep("min",num.factors))
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   edge.label <- rep("",num.var*k)
   edge.name <- rep("",num.var*k)
   names(edge.label) <-  seq(1:num.var*k) 
   edge.dir <- rep("forward",num.var*k)        
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
       
 if(!is.null(Phi)) {
    k <- num.var +1 
    
    for (f in 2:num.factors)  {
     for (f1 in 1:(f-1)) { if(Phi[f,f1] > cut) {
        clust.graph <- addEdge(fact[f1], fact[f], clust.graph,1) 
        edge.label[k] <- round(Phi[f,f1],digits)
        edge.name[k] <- paste(fact[f1],"~",fact[f],sep="")
        edge.dir[k] <- paste("both")
        k <- k+1}
                         }
                               } 
    }
                  
 nAttrs <- list()  #node attributes
 eAttrs <- list()  #edge attributes


 if (!is.null(labels)) {var.labels <- c(labels,fact)
  names(var.labels) <-  nodes(clust.graph)
  nAttrs$label <- var.labels
  names(edge.label) <- edge.name
  } 
    names(edge.label) <- edge.name
    names(edge.dir) <- edge.name
 
 nAttrs$shape <- graph.shape
 nAttrs$rank <- graph.rank
 eAttrs$label <- edge.label
 eAttrs$dir <- edge.dir
 #eAttrs$font <- edge.font
 attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE),graph=list(rankdir=rank.direction, fontsize=edge.font[2],bgcolor="white" ))
 obs.var <- subGraph(vars,clust.graph)
 cluster.vars <- subGraph(fact,clust.graph)
 observed <- list(list(graph=obs.var,cluster=TRUE,attrs=c(rank="sink")),list(graph=cluster.vars,cluster=FALSE ,attrs=c(rank = "source"))) #this crashes for correlated factors solution

 observed <- list(list(graph=obs.var,cluster=TRUE,attrs=c(rank="sink")))   #this does not lead to a crash
plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed,main=main)  

 if(!is.null(out.file) ){toDotty(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs) }


return(clust.graph)
   }
 
