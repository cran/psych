#factor analysis and sem diagrams
#based upon fa.graph with some ideas taken from the diagram and shape packages of Karline Soetaert 
#version of September 30, 2009
#developed to replace Rgraphviz which is too much of a pain to install
#Rgraphviz uses a NEL  (Node, Edge, Label) represenation while diagram uses a complete linking matrix
#thus, I am trying to combine these two approaches

"fa.diagram" <-
  function(fa.results,sort=TRUE,labels=NULL,cut=.3,simple=TRUE,errors=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main="Factor Analysis",cex=NULL, ...) {
   col <- c("black","red")
   if(is.null(cex)) cex <- 1
  Phi <- NULL  #the default case
 if(sort) fa.results <- fa.sort(fa.results) 
 if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {factors <- as.matrix(fa.results$loadings)
                if(!is.null(fa.results$Phi)) Phi <- fa.results$Phi} else {factors <- fa.results}
   
       nvar <- dim(factors)[1]   #how many variables?
   if (is.null(nvar) ){nvar <- length(factors)
       num.factors <- 1} else {
         num.factors <- dim(factors)[2]}
#first some basic setup parameters 
  
   nvar <- dim(factors)[1]   #how many variables?
   if (is.null(nvar) ){nvar <- length(factors)
       num.factors <- 1} else {
         num.factors <- dim(factors)[2]}
   
   if (is.null(rownames(factors))) {rownames(factors) <- paste("V",1:nvar,sep="") }
   if (is.null(colnames(factors))) {colnames(factors) <- paste("F",1:num.factors,sep="") }
   
   var.rect <- list()
   fact.rect <- list()
   max.len <- max(nchar(rownames(factors)))*rsize
  plot(0,type="n",xlim=c(-max.len/2,nvar+1),ylim=c(0,nvar+1),asp=1,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
  cex <-  min(cex,20/nvar)
  for (v in 1:nvar) { 
 	 var.rect[[v]] <- dia.rect(1,nvar-v+1,rownames(factors)[v],xlim=c(0,nvar),ylim=c(0,nvar),cex=cex,...)
     }
   f.scale <- (nvar+ 1)/(num.factors+1)
   f.shift <- nvar/num.factors
   for (f in 1:num.factors) {
   		fact.rect[[f]] <- dia.ellipse(nvar-2,(num.factors+1-f)*f.scale,colnames(factors)[f],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
     		for (v in 1:nvar)  {
     		
    			if(simple && (abs(factors[v,f]) == max(abs(factors[v,])) )  && (abs(factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) { 
    			dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =round(factors[v,f],digits),col=((sign(factors[v,f])<0) +1))
    			

    }
   }
   }
   
   if(!is.null(fa.results$Phi)) { for (i in 2:num.factors) {
     for (j in 1:(i-1)) {
     if(Phi[i,j] > cut) {
        dia.curve(from=c(nvar-2+ e.size*nvar,(num.factors+1-i)*f.scale),to=c(nvar -2+ e.size*nvar,(num.factors+1-j)*f.scale),labels=round(Phi[i,j],digits),scale=(i-j),...)}
						}
  															 }
 
						}
	
  if (errors) {for (v in 1:nvar) {
       dia.self(location=var.rect[[v]],scale=.5,side=side)  }
       }
}                  