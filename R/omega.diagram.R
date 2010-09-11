#factor analysis and sem diagrams
#based upon fa.graph with some ideas taken from the diagram and shape packages of Karline Soetaert 
#version of September 20, 2009
#developed to replace Rgraphviz which is too much of a pain to install
#Rgraphviz uses a NEL  (Node, Edge) represenation while diagram uses a complete linking matrix
#thus, I am trying to combine these two approaches

"omega.diagram" <-
  function(om.results,sl=TRUE,sort=TRUE,labels=NULL,cut=.2,gcut=.2,simple=TRUE,errors=FALSE,
    digits=1,e.size=.1,rsize=.15,side=3,main=NULL,cex=NULL,color.lines=TRUE, ...) {
     if(color.lines) { colors <- c("black","red")} else {colors <- c("black","black") }
  Phi <- NULL  #the default case
  if(is.null(cex)) cex <- 1
  
 
 if((length(class(om.results)) > 1)  && ( (class(om.results)[2] == "omegaSem")) ) {
      #  factors <- as.matrix(om.results$omega.efa$cfa.loads[,2:ncol(om.results$omega.efa$cfa.loads)])
      if(sort) om.results$omega.efa$cfa.loads <- fa.sort(om.results$omega.efa$cfa.loads)
      factors <- as.matrix(om.results$omega.efa$cfa.loads)
        gloading <- om.results$omega.efa$cfa.loads[,1,drop=FALSE]
        nvar <- num.var <- nrow(gloading)
        num.factors <- ncol(factors) -1
        sl=TRUE
        main <- "Omega from SEM" 
    
        } else {
 if(sort) om.results <- fa.sort(om.results)   #usually sort, but sometimes it is better not to do so
 if (sl) {factors <- as.matrix(om.results$schmid$sl) 
         if(is.null(main)) {main <- "Omega with Schmid Leiman Transformation" }
         } else {factors <- as.matrix(om.results$schmid$oblique)
          if(is.null(main)) {main <- "Hierarchical (multilevel) Structure" }
         }
       gloading <- om.results$schmid$gloading
        nvar <- num.var <- dim(factors)[1]   #how many variables?
   if (sl) {num.factors <- dim(factors)[2] -4 } else {num.factors <- dim(factors)[2]
       }
 }
   
   e.size <- e.size * 10/ nvar   #this is an arbitrary setting that seems to work
#first some basic setup parameters 
  
   
   
    vars <- paste("V",1:num.var,sep="")  
   if (!is.null(labels)) {vars <- paste(labels)} else{vars <- rownames(factors) }
   if(sl) {fact <- c("g",paste("F",1:num.factors,"*",sep="")) } else {fact <- c("g",paste("F",1:num.factors,sep="")) }   # e.g.  "g"  "F'1" "F2" "F3"
   var.rect <- list()
   fact.rect <- list()
    max.len <- max(nchar(rownames(factors)))*rsize
   cex <-  min(cex,20/nvar)
  plot(0,type="n",xlim=c(-max.len/2,nvar+1),ylim=c(1,nvar+1),asp=1,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
   if (sl) {vloc <- (nvar+1)/2
            gloc <- 1
            grouploc <- nvar
            start <- 1
            end <- num.factors+1} else {
            vloc <- 1
            gloc <- nvar
            grouploc <- (nvar+1)/2
            start <- 0
            end <- num.factors
            }
  for (v in 1:nvar) { 
 	 var.rect[[v]] <- dia.rect(vloc,nvar-v+1,rownames(factors)[v],xlim=c(0,nvar),ylim=c(0,nvar),cex=cex,...)
     }
   f.scale <- (nvar+ 1)/(num.factors+1)
   f.shift <- nvar/num.factors
   for (f in 1:num.factors) {
   		fact.rect[[f]] <- dia.ellipse(grouploc,(num.factors+1-f)*f.scale,colnames(factors)[f+start],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
     		for (v in 1:nvar)  {
    			if (abs(factors[v,f+start]) > cut) {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,col=colors[((sign(factors[v,f+start])<0) +1)],lty=((sign(factors[v,f+start])<0)+1),labels=round(factors[v,f+start],digits))
                               }
                             }
                     }
   
  g.ellipse <-  dia.ellipse(gloc,(num.var+1)/2,"g",xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
   if(!sl) { 
   for (f in 1:num.factors) {
              dia.arrow(from=g.ellipse,to=fact.rect[[f]],col=colors[((sign(gloading[f])<0) +1)],lty=((sign(gloading[f])<0) +1),labels=round(gloading[f],digits))                
                            }
             } else {
              for (i in 1:nvar) {
              if(abs(factors[i,1]) > gcut) {
               dia.arrow(from=g.ellipse,to=var.rect[[i]]$left,col=colors[((sign(factors[i,1])<0) +1)],lty=((sign(factors[i,1])<0)+1),labels=round(factors[i,1],digits))}
                                 }
 
	  }
  if (errors) {for (v in 1:nvar) {
       dia.self(location=var.rect[[v]],scale=.5,side=side)  }
       }
}                  