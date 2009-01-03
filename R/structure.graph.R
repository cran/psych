#Created April 4, 2008 

#in process  January 2, 2009

"structure.graph" <-
function(xmodel,ymodel=NULL,Phi=NULL, out.file=NULL,labels=NULL,cut=.3,simple=TRUE,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10),  rank.direction=c("RL","TB","LR","BT"), digits=1,title="Structural model", ...){
    	if (!require(Rgraphviz)) {stop("I am sorry, you need to have loaded the Rgraphviz package")}
 if(!is.matrix(xmodel) & !is.data.frame(xmodel))  {
        if(!is.null(xmodel$Phi)) Phi <- xmodel$Phi
        xmodel <- as.matrix(xmodel$loadings)
      } else {xmodel <- xmodel}   
 if(!is.matrix(xmodel) ) {factors <- as.matrix(xmodel)} else {factors <- xmodel}
 
 
 
   rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
  num.y <- 0   #we assume there is nothing there
  num.var <- dim(factors)[1]   #how many x variables?
  if (is.null(num.var) ){num.var <- length(factors)
                          num.factors <- 1} else {
            num.factors <- dim(factors)[2]}
   
   num.xfactors <- num.factors
   num.xvar <- num.var
   vars <- rownames(xmodel)
  if(is.null(vars) ) {vars <- paste("Xv",1:num.var,sep="")  }
  fact <- colnames(xmodel)
  if (is.null(fact)) { fact <- paste("X",1:num.factors,sep="") }
  
   if (!is.null(ymodel)) { 
   if(!is.matrix(ymodel) & !is.data.frame(ymodel))  {ymodel <- as.matrix(ymodel$loadings)} else {ymodel <- ymodel}   
 if(!is.matrix(ymodel) ) {y.factors <- as.matrix(ymodel)} else {y.factors <- ymodel}
   
   num.y <- dim(y.factors)[1] 
      if (is.null(num.y)) {
         num.y <- length(ymodel)
         num.yfactors <- 1} else {
         num.yfactors <- dim(ymodel)[2]
        }  
      
     yvars <- rownames(ymodel)
     if(is.null(yvars)) {yvars <- paste("Yv",1:num.y,sep="")  }
     vars <- c(vars,yvars)
     
      yfact <- colnames(ymodel)
      if(is.null(yfact)) {yfact <- paste("Y",1:num.yfactors,sep="") }
      fact <- c(fact,yfact)
     num.var <- num.xvar + num.y
     num.factors <- num.xfactors + num.yfactors
     }
     
   
   #now draw the x part
   k <- num.factors  
      
   clust.graph <-  new("graphNEL",nodes=c(vars,fact),edgemode="directed")
   graph.shape <- c(rep("box",num.var),rep("ellipse",num.factors))  #define the shapes
   if (num.y > 0) {graph.rank <- c(rep("min",num.xvar),rep("max",num.y),rep("",num.factors))} else {
                  graph.rank <- c(rep("min",num.var),rep("",num.factors))} 
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   edge.label <- rep("",num.var*k)  #makes too many, but in case of a fully saturated model, this might be necessary
   edge.name <- rep("",num.var*k)
   names(edge.label) <-  seq(1:num.var*k) 
   edge.dir <-rep("forward",num.var*k)
   names(edge.dir) <-rep("",num.var*k)
   edge.arrows <-rep("open",num.var*k)
   names(edge.arrows) <-rep("",num.var*k)
  #show the cluster structure with ellipses

  l <- factors   # ?
  if (num.xfactors ==1) { 
       
    for (i in 1:num.xvar) { clust.graph <- addEdge(fact[1], vars[i], clust.graph,1) 
                                           edge.label[i] <- round(factors[i],digits)
                                           edge.name[i] <- paste(fact[1],"~",vars[i],sep="")
                        }  
                        k <- num.var+1 
      } else { 
        #all loadings > cut in absolute value
                   k <- 1
                   for (i in 1:num.xvar) {
                   for (f in 1:num.xfactors) { if (abs(factors[i,f]) > cut) {
               clust.graph <- addEdge(fact[f], vars[i], clust.graph,1) 
                              edge.label[k] <- round(factors[i,f],digits)
                              edge.name[k] <- paste(fact[f],"~",vars[i],sep="")
                         k <- k+1 }   #end of if 
                         }  
                        }      
       }   #end of if num.xfactors ==1 
       
  #now, if there is a ymodel, do it for y model 
  
  if(!is.null(ymodel)) { 
  if (num.yfactors ==1) {     
    for (i in 1:num.y) { clust.graph <- addEdge( yvars[i],fact[1+num.xfactors], clust.graph,1) 
                         edge.label[k] <- round(y.factors[i],digits)
                         edge.name[k] <- paste(yvars[i],"~",fact[1+num.xfactors],sep="")
                         edge.dir[k] <- paste("back")
                         k <- k +1
                        }                      
      } else { 
        #all loadings > cut in absolute value
                   for (i in 1:num.y) {
                   for (f in 1:num.yfactors) { if (abs(y.factors[i,f]) > cut) {clust.graph <- addEdge( vars[i+num.xvar],fact[f+num.xfactors], clust.graph,1) 
                         edge.label[k] <- round(y.factors[i,f],digits)
                         edge.name[k] <- paste(vars[i+num.xvar],"~",fact[f+num.xfactors],sep="")
                          edge.dir[k] <- paste("back")
                         k <- k+1 }   #end of if 
                         }  #end of factor
                        } # end of variable loop
           
       }   #end of if num.yfactors ==1 
  
  }   #end of if.null(ymodel)
                
 nAttrs <- list()  #node attributes
 eAttrs <- list()  #edge attributes


 if (!is.null(labels)) {var.labels <- c(labels,fact)
  names(var.labels) <-  nodes(clust.graph)
  nAttrs$label <- var.labels
  names(edge.label) <- edge.name
  } 
  
  if(!is.null(Phi)) {for (i in 2:num.factors) {
                         for (j in 1:(i-1)) { if(abs(Phi[i,j]) > cut ) {
                            clust.graph <- addEdge( fact[i],fact[j],clust.graph,1)
                            edge.label[k] <- round(Phi[i,j],digits)
                            edge.name[k] <- paste(fact[i],"~",fact[j],sep="")
                            if(Phi[i,j] != Phi[j,i]){edge.dir[k] <- "back"} else {edge.dir[k] <- "both"}}
                          k <- k + 1
                         }
                         }
                    
   
  
  
  }
  
  obs.xvar <- subGraph(vars[1:num.xvar],clust.graph)
if (!is.null(ymodel)) {obs.yvar <- subGraph(vars[(num.xvar+1):num.var],clust.graph)}
#obs.var <- subGraph(vars,clust.graph)
 cluster.vars <- subGraph(fact,clust.graph)
 #observed <- list(list(graph=obs.xvar,cluster=TRUE,attrs=c(rank="min")))
 
 
 names(edge.label) <- edge.name
 names(edge.dir) <- edge.name
 names(edge.arrows) <- edge.name
 nAttrs$shape <- graph.shape
 nAttrs$rank <- graph.rank
 eAttrs$label <- edge.label
 eAttrs$dir<- edge.dir
  eAttrs$arrowhead <- edge.arrows
 eAttrs$arrowtail<- edge.arrows
 attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE),graph=list(rankdir=rank.direction, fontsize=10,bgcolor="white" ))
 
 
if (!is.null(ymodel)) {observed <- list(list(graph=obs.xvar,cluster=TRUE,attrs=c(rank="sink")),list(graph=obs.yvar,cluster=TRUE,attrs=c(rank="source"))) } else {
                       observed <- list(list(graph=obs.xvar,cluster=TRUE,attrs=c(rank="max")))}
 plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed,main=title) 
if(!is.null(out.file) ){toDot(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs) }

#return(list(nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs))  #useful for debugging 

   }
   
 
