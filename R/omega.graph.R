#Created May 20, 2007
#modified June 2 to clarify the Rgraphviz issue
#modified July 12 to fix label problem
#take the output from omega and graph it
"omega.graph" <-
function(om.results,out.file=NULL,sl=TRUE,labels=NULL,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10), rank.direction="RL", digits=1,title="Omega", ...){
   
     if(!require(Rgraphviz)) {stop("I am sorry, you need to have the  Rgraphviz package installed")}
    
   if (sl) {factors <- as.matrix(om.results$schmid$sl)   } else{factors <- as.matrix(om.results$schmid$oblique)}
   rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
  
   num.var <- dim(factors)[1]   #how many variables?
  if (sl) {num.factors <- dim(factors)[2] -3 } else {num.factors <- dim(factors)[2]}
   gloading <- om.results$schmid$gloading
   vars <- paste("V",1:num.var,sep="")  
   if (!is.null(labels)) {vars <- paste(labels,1:num.var,sep=" ")}
   fact <- c("g",paste("F",1:num.factors,sep=""))
   clust.graph <-  new("graphNEL",nodes=c(vars,fact),edgemode="directed")
   graph.shape <- c(rep("box",num.var),rep("ellipse",num.factors+1))
   graph.rank <- c(rep("sink",num.var),rep("",num.factors+1))
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   if (sl) {edge.label <- rep("",num.var*2)
            edge.name <- rep("",num.var*2)
            names(edge.label) <-  seq(1:num.var*2) } else {
             edge.label <- rep("",num.var+num.factors)
            edge.name <- rep("",num.var+num.factors)
            names(edge.label) <-  seq(1:num.var+num.factors) }
            
  #show the cluster structure with ellipses
   if (sl) {
   l <- matrix(factors[,2:(num.factors+1)],ncol=num.factors) } else { l <- factors }
   m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max), 
        ncol = 1)
        
  if (sl) { k <- num.var
         for (i in 1:num.var) {  clust.graph <- addEdge(fact[1], vars[i], clust.graph,1) 
             edge.label[i] <- round(factors[i,1],digits)
              edge.name[i] <- paste(fact[1],"~",vars[i],sep="")
             } 
         } else { k <- num.factors+1
          for (j in 1:num.factors) {clust.graph <- addEdge(fact[1], fact[j+1], clust.graph,1) 
          edge.label[j] <- round(gloading[j],digits)
          edge.name[j] <- paste(fact[1],"~",fact[j+1],sep="")
         } 
                 }
   for (i in 1:num.var) {  clust.graph <- addEdge(fact[1+m1[i]], vars[i], clust.graph,1) 
                         edge.label[i+k] <- round(l[i,m1[i]],digits)
                         edge.name[i+k] <- paste(fact[1+m1[i]],"~",vars[i],sep="")
                        }  

 if(FALSE) {      
     edge.label[(i-1)*2+1] <- results[i,"r1"]
     edge.name [(i-1)*2+1]  <- paste(row.names(results)[i],"~", results[i,1],sep="")
      
     } 
     
 nAttrs <- list()  #node attributes
 eAttrs <- list()  #edge attributes

if(FALSE) {
 if (!is.null(labels)) {var.labels <- c(labels,fact)
  							names(var.labels) <-  nodes(clust.graph)
 							nAttrs$label <- var.labels
  							names(edge.label) <- edge.name
  						} 
  					}
				
 names(edge.label) <- edge.name
 nAttrs$shape <- graph.shape
 nAttrs$rank <- graph.rank
 eAttrs$label <- edge.label
 attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE),graph=list(rankdir="RL", fontsize=10,bgcolor="white" ))
 obs.var <- subGraph(vars,clust.graph)
 cluster.vars <- subGraph(fact,clust.graph)
 observed <- list(list(graph=obs.var,cluster=TRUE,attrs=c(rank="")))
 plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed,main=title) 
if(!is.null(out.file) ){toDot(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs) }
return(clust.graph)
   }
 
