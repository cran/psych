#modified January 20, 2009 to create sem commands
#modified May 30, 2008 to try to get arrows going both ways in the sl option.
#this now works, but sometimes two lines don't have arrows.
#Created May 20, 2007
#modified June 2 to clarify the Rgraphviz issue
#modified July 12 to fix label problem
#take the output from omega and graph it
#fixed Sept 16 to draw sl solutions correctly
"omega.graph" <-
function(om.results,out.file=NULL,sl=TRUE,labels=NULL,
   size=c(8,6), node.font=c("Helvetica", 14),
    edge.font=c("Helvetica", 10),  rank.direction=c("RL","TB","LR","BT"), digits=1,title="Omega", ...){
   
     
     if(!requireNamespace('Rgraphviz')) {stop("I am sorry, you need to have the  Rgraphviz package installed")
        #create several dummy functions to get around the "no visible global function definition" problem
    	nodes <- function() {}
    	addEdge <- function() {}
    	subGraph <- function(){} }
    # if(!requireNamespace(graph)) {stop("I am sorry, you need to have the  graph package installed") }
    
   if (sl) {factors <- as.matrix(om.results$schmid$sl)   } else{factors <- as.matrix(om.results$schmid$oblique)}
   rank.direction <- match.arg(rank.direction)
  #first some basic setup parameters 
  
   num.var <- dim(factors)[1]   #how many variables?
  if (sl) {num.factors <- dim(factors)[2] -4 } else {num.factors <- dim(factors)[2]}
   gloading <- om.results$schmid$gloading
   vars <- paste("V",1:num.var,sep="")  
   if (!is.null(labels)) {vars <- paste(labels)} else{vars <- rownames(factors) }
   if(sl) {fact <- c("g",paste("F",1:num.factors,"*",sep="")) } else {fact <- c("g",paste("F",1:num.factors,sep="")) }   # e.g.  "g"  "F'1" "F2" "F3"
  clust.graph <-  new("graphNEL",nodes=c(vars,fact),edgemode="directed")
   graph.shape <- c(rep("box",num.var),rep("ellipse",num.factors+1))
   graph.rank <- c("sink", rep("same",num.var),rep("source",num.factors)) #this doesn't seem to make a difference
   names(graph.shape) <- nodes(clust.graph)
   names(graph.rank) <- nodes(clust.graph)
   
   if (sl) {edge.label <- rep("",num.var*2)    #this basically just sets up the vectors to be the right size
            edge.dir <-rep("forward",num.var*2)
           # edge.arrows <-rep("open",num.var+num.factors)
           edge.arrows <-rep("open",num.var*2)
            edge.name <- rep("",num.var*2)
            names(edge.label) <-  seq(1:num.var*2)
            names(edge.dir) <-rep("",num.var*2)
            #names(edge.arrows) <-rep("",num.var+num.factors)
            names(edge.arrows) <-rep("",num.var*2)
            sem <- matrix(rep(NA,6*(2*num.var + num.factors)),ncol=3)  #used for sem diagram
            } else {
            
            edge.label <- rep("",num.var+num.factors)
            edge.name <- rep("",num.var+num.factors)
            edge.arrows <-rep("open",num.var+num.factors)
            edge.dir <-rep("forward",num.var*2)
            names(edge.label) <-  seq(1:num.var+num.factors)
            names(edge.dir) <-  seq(1:num.var+num.factors)
            names(edge.arrows) <-  seq(1:num.var+num.factors)
            sem <- matrix(rep(NA,6*(num.var + num.factors)+3),ncol=3)  #used for sem diagram
                       }
            
  #show the cluster structure with ellipses
   if (sl) {
   l <- matrix(factors[,2:(num.factors+1)],ncol=num.factors) } else { l <- factors }
   
   m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max),  ncol = 1)
        
  if (sl) { 
          k <- num.var
         for (i in 1:num.var) {  clust.graph <- addEdge( vars[i],fact[1], clust.graph,1) 
                                 edge.label[i] <- round(factors[i,1],digits)
                                 edge.name[i] <- paste(vars[i],"~",fact[1],sep="")
                                 edge.arrows[i] <- paste("open")
                                 edge.dir[i] <- paste("back")
                                 sem[i,1] <- paste(fact[1],"->",vars[i],sep="")
                                 sem[i,2] <- vars[i]               
            } 
            
         } else { 
         
          k <- num.factors
          for (j in 1:num.factors) {clust.graph <- addEdge(fact[1], fact[j+1], clust.graph,1)  #hierarchical g 
          edge.label[j] <- round(gloading[j],digits)
          edge.name[j] <- paste(fact[1],"~",fact[j+1],sep="")
           sem[j,1] <- paste(fact[1],"->",fact[1+j],sep="")
           sem[j,2] <- paste("g",fact[1+j],sep="")
                 } 
                 }
   for (i in 1:num.var) {  clust.graph <- addEdge(fact[1+m1[i]], vars[i], clust.graph,1) 
                         edge.label[i+k] <- round(l[i,m1[i]],digits)
                         edge.name[i+k] <- paste(fact[1+m1[i]],"~",vars[i],sep="")
                         edge.arrows[i+k] <- paste("open")
                           sem[i+k,1] <- paste(fact[1+m1[i]],"->",vars[i],sep="")
                          sem[i+k,2] <- paste(fact[1+m1[i]],vars[i],sep="")
                        }  

     
#    edge.label[(i-1)*2+1] <- results[i,"r1"]
#    edge.name [(i-1)*2+1]  <- paste(row.names(results)[i],"~", results[i,1],sep="")
 if(sl) {
       k <- num.var*2
      for (i in 1:num.var) {
            sem[i+k,1] <- paste(vars[i],"<->",vars[i],sep="")
            sem[i+k,2] <- paste("e",i,sep="")
                            }
        k <- k + num.var
       for (f in 1:num.factors) {
             sem[f+k,1] <- paste(fact[1+f],"<->",fact[1+f],sep="")
             sem[f+k,3] <- "1"
                                 }
        k <- k+ num.factors
           sem[k+1,1] <- paste("g <->g")
           sem[k+1,3] <- "1"
           k<- k+1
           
          } else { 
          
          k <- num.var + num.factors
          for (i in 1:num.var) {
            sem[i+k,1] <- paste(vars[i],"<->",vars[i],sep="")
            sem[i+k,2] <- paste("e",i,sep="")
                            }
            k <- 2*num.var + num.factors
          for (f in 1:num.factors) {
            sem[f+k,1] <- paste(fact[f+1],"<->",fact[f+1],sep="")
            sem[f+k,3] <- "1"
                            }
             k <- 2*num.var + 2*num.factors
             sem[k+1,1] <- paste("g<->g")
             sem[k+1,3] <- "1"
             k <- k+1
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
 names(edge.dir) <- edge.name
 names(edge.arrows) <- edge.name
 nAttrs$shape <- graph.shape
 nAttrs$rank <- graph.rank
 eAttrs$label <- edge.label
if(sl) { eAttrs$dir<- edge.dir
         
         eAttrs$arrowhead <- edge.arrows
         eAttrs$arrowtail<- edge.arrows
       
  } 

 attrs <- list(node = list(shape = "ellipse", fixedsize = FALSE),graph=list(rankdir=rank.direction, fontsize=10,bgcolor="white" ))
 
# obs.var <- subGraph(vars,clust.graph)
# cluster.vars <- subGraph(fact,clust.graph)
# observed <- list(list(graph=obs.var,cluster=TRUE,attrs=c(rank="")))
# plot(clust.graph, nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs,subGList=observed,main=title) 


plot(clust.graph,nodeAttrs = nAttrs,edgeAttrs = eAttrs,attrs = attrs,main=title)   #not clear if the subGList makes any difference
if(!is.null(out.file) ){toDotty(clust.graph,out.file,nodeAttrs = nAttrs, edgeAttrs = eAttrs, attrs = attrs) }
colnames(sem) <- c("Path","Parameter","Initial Value")
return(sem=sem[1:k,])
   }
 
