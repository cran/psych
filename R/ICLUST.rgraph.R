#create Rgraphviz commands
#create dot file for graphviz from ICLUST output
#modified from ICLUST.graph to produce Rgraphviz output

"ICLUST.rgraph"  <- 
function(ic.results,out.file = NULL, min.size=1,short=FALSE,labels=NULL,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction="RL", digits=2,title="ICLUST", ...){
 if(!require(Rgraphviz)) {stop("I am sorry, you need to have the  Rgraphviz package installed")}
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
 plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed,main=title) 
 if (!is.null(out.file)) {toDot(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed) }
   }
 
#test with Harman data sets
#ic1 <- ICLUST(Harman74.cor$cov)
#ic4 <- ICLUST(Harman74.cor$cov,nclusters=4)
#ic8 <- ICLUST(Harman23.cor$cov)
#ICLUST.rgraph(ic1)
#ICLUST.rgraph(ic4)
#ICLUST.rgraph(ic8,labels=colnames(Harman23.cor$cov))
 