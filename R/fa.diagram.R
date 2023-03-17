#factor analysis and sem diagrams
#based upon fa.graph with some ideas taken from the diagram and shape packages of Karline Soetaert 
#version of September 30, 2009
#developed to replace Rgraphviz which is too much of a pain to install
#Rgraphviz uses a NEL  (Node, Edge, Label) representation while diagram uses a complete linking matrix
#thus, I am trying to combine these two approaches
#some small revisions, March 2015 to allow cex to be passed to dia.ellipse and dia.rect
#revised June, 2020 to use dia.ellipse and dia.rect and dia.arrow in vector mode


"fa.diagram" <-
  function(fa.results,Phi=NULL,fe.results=NULL,sort=TRUE,labels=NULL,cut=.3,simple=TRUE,errors=FALSE,g=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main,cex=NULL,l.cex=NULL,marg=c(.5,.5,1,.5),adj=1,ic=FALSE, ...) {
    if(length(class(fa.results)) > 1) {if(inherits(fa.results, 'principal')) {pc <- TRUE} else {pc <- FALSE}} else { pc <- FALSE}
    if(ic) pc <- TRUE
    old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back
   col <- c("black","red")
   if(missing(main)) if(is.null(fe.results)) {if(pc) {main <- "Components Analysis" } else {main <- "Factor Analysis" }} else {main <- "Factor analysis and extension"}
   if(!is.matrix(fa.results) && !is.null(fa.results$fa) && is.list(fa.results$fa)) fa.results <- fa.results$fa
   if(is.null(cex)) cex <- 1
   if(is.null(l.cex)) l.cex <- 1
  #Phi <- NULL  #the default case
 if(sort) { if(g) {temp  <- fa.sort(fa.results[,-1])
                temp2 <- fa.results[,1]
          fa.results <- cbind(g=temp2[rownames(temp)],temp)
           } else {fa.results <- fa.sort(fa.results)}    #if we have g loadings, don't sort the entire array
          if(!is.null(fe.results)) { fe.results <- fa.sort(fe.results)} }
 if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {factors <- as.matrix(fa.results$loadings)
                if(!is.null(fa.results$Phi)) {Phi <- fa.results$Phi} else {
                       if(!is.null(fa.results$cor)) {Phi<- fa.results$cor} 
                       }} else {factors <- fa.results}
   
       nvar <- dim(factors)[1]   #how many variables?
   if (is.null(nvar) ){nvar <- length(factors)
       num.factors <- 1} else {
         num.factors <- dim(factors)[2]}
#first some basic setup parameters 
  
   nvar <- dim(factors)[1]   #how many variables?
   e.size = e.size*16*cex/nvar
   if (is.null(nvar) ){nvar <- length(factors)
       num.factors <- 1} else {
         num.factors <- dim(factors)[2]}
   
   if (is.null(rownames(factors))) {rownames(factors) <- paste("V",1:nvar,sep="") }
   if (is.null(colnames(factors))) {colnames(factors) <- paste("F",1:num.factors,sep="") }
   
   var.rect <- list()
   fact.rect <- list()

   
   max.len <- max(nchar(rownames(factors)))*rsize
   x.max <-  max((nvar+1),6) 
  
   limx=c(-max.len/2,x.max)
   n.evar <- 0

    if(!is.null(fe.results)) {n.evar <- dim(fe.results$loadings)[1]
    limy <- c(0,max(nvar+1,n.evar+1))} else {
     limy=c(0,nvar+1) }
     top <- max(nvar,n.evar) + 1
# plot(0,type="n",xlim=limx,ylim=limy,asp=1,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main,...)
  plot(0,type="n",xlim=limx,ylim=limy,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main,...)
 
   max.len <- max(strwidth(rownames(factors)),strwidth("abc"))/1.8  #slightly more accurate, but needs to be called after plot is opened
    limx=c(-max.len/2,x.max)  
     cex <-  min(cex,20/x.max)
 if(g) {left <- .3*x.max     #where should the variable boxes go?  It depends upon g
        middle <- .6*x.max
        gf <- 2 } else {left <- 0
        middle <- .5*x.max
        gf <- 1}  
 for (v in 1:nvar) { 
 	 var.rect[[v]] <- dia.rect(left,top -v - max(0,n.evar-nvar)/2  ,rownames(factors)[v],xlim=limx,ylim=limy,cex=cex,draw=FALSE,...)  #find their locations but don't draw them
     }
     all.rects.x <- rep(left,nvar)
     all.rects.y <-  top - 1:nvar -  max(0,n.evar-nvar)/2
     all.rects.names <- rownames(factors)[1:nvar]
    dia.rect(all.rects.x, all.rects.y,all.rects.names)    #now draw them as a vector
    
   f.scale <- (top)/(num.factors+1)
   f.shift <- max(nvar,n.evar)/num.factors
   if(g) {fact.rect[[1]] <- dia.ellipse(-max.len/2,top/2,colnames(factors)[1],xlim=limx,ylim=limy,e.size=e.size,cex=cex,...)
          	for (v in 1:nvar)  {if(simple && (abs(factors[v,1]) == max(abs(factors[v,])) )  && (abs(factors[v,1]) > cut) | (!simple && (abs(factors[v,1]) > cut))) { 
    			dia.arrow(from=fact.rect[[1]],to=var.rect[[v]]$left,labels =round(factors[v,1],digits),col=((sign(factors[v,1])<0) +1),lty=((sign(factors[v,1])<0)+1))
    	 }}}
    	  
    	  text.values <- list()
    	  tv.index <- 1
   for (f in gf:num.factors) {  #body  34
   		if (pc) {fact.rect[[f]] <- dia.rect(left+middle,(num.factors+gf-f)*f.scale,colnames(factors)[f],xlim=limx,ylim=limy,cex=cex,draw=FALSE,...) 
   		} else {fact.rect[[f]] <- dia.ellipse(left+middle,(num.factors+gf-f)*f.scale,colnames(factors)[f],xlim=limx,ylim=limy,e.size=e.size,cex=cex,draw=FALSE,...)}
     		
     		
     		for (v in 1:nvar)  {
     		
    			if(simple && (abs(factors[v,f]) == max(abs(factors[v,])) )  && (abs(factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) { 
    		if(pc) {text.values[[tv.index]] <- dia.arrow(to=fact.rect[[f]],from=var.rect[[v]]$right,labels =round(factors[v,f],digits),col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0)+1),adj=f %% adj ,cex=cex,draw=FALSE)
    		    tv.index <- tv.index + 1 
    		} else {text.values[[tv.index]] <- dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =round(factors[v,f],digits),col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0)+1),adj=f %% adj +1,cex=cex,draw=FALSE)
    		       tv.index <- tv.index + 1  }   
    			 }
   }
   }
    #now draw all the factors/components and their associated loadings
    #first, just the factors and their labels
    
     
      tv <- matrix(unlist(fact.rect),nrow=num.factors,byrow=TRUE)   #unfortunately, this loses whether it is an elipse or not
      all.rects.x <- tv[,5] #the center of the figure
      all.rects.y <- tv[,2]
      all.rects.names <- colnames(factors)
     if(pc) {
      dia.rect(all.rects.x, all.rects.y,all.rects.names)} else {dia.multi.ellipse(all.rects.x,all.rects.y,all.rects.names,e.size=e.size)
      }  #we need to remember whether to draw ellipses or boxes
     
   #now show all the loadings
   
   tv <- matrix(unlist(text.values),byrow=TRUE,ncol=21)
        #text(tv[,1],tv[,2],tv[,3],cex=tv[,5])
        text(tv[,1],tv[,2],tv[,3],cex=l.cex)
        arrows(x0=tv[,6],y0=tv[,7],x1=tv[,8],y1=tv[,9],length=tv[1,10],angle=tv[1,11],code=1,col=tv[,20],lty=tv[,21])
        arrows(x0=tv[,13],y0=tv[,14],x1=tv[,15],y1=tv[,16],length=tv[1,17],angle=tv[1,18],code=2,col=tv[,20],lty=tv[,21])



   if(!is.null(Phi) && (ncol(Phi) >1)) { 
       curve.list <- list()
       for (i in 2:num.factors) {
  
     for (j in 1:(i-1)) {
     if(abs(Phi[i,j]) > cut) {
       # dia.curve(from=c(x.max-2+ e.size*nvar,(num.factors+1-i)*f.scale),to=c(x.max -2+ e.size*nvar,(num.factors+1-j)*f.scale),labels=round(Phi[i,j],digits),scale=(i-j),...)}
		d.curve <- dia.curved.arrow(from=fact.rect[[j]]$right,to=fact.rect[[i]]$right,labels=round(Phi[i,j],digits),scale=(i-j),draw=FALSE,cex=cex,l.cex=l.cex,...)
		curve.list <- c(curve.list,d.curve)
		} }
  															 }
  			multi.curved.arrow(curve.list,l.cex,...)												 
 
						}
	self.list <- list()
  if (errors) {for (v in 1:nvar) {
      d.self <-  dia.self(location=var.rect[[v]],scale=.5,side=side)
      self.list <- c(self.list,d.self)  }
       }
    if(length(self.list) > 0)  multi.self(self.list)
       
   if(!is.null(fe.results)) {
     e.loadings <- fe.results$loadings
     
     for (v in 1:n.evar) { 
 	 var.rect[[v]] <- dia.rect(x.max,top-v-max(0,nvar-n.evar)/2,rownames(e.loadings)[v],xlim=limx,ylim=limy,cex=cex,...)
 	 for(f in 1:num.factors) {
 	 if(simple && (abs(e.loadings[v,f]) == max(abs(e.loadings[v,])) )  && (abs(e.loadings[v,f]) > cut) | (!simple && (abs(e.loadings[v,f]) > cut))) { 
    			dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$left,labels =round(e.loadings[v,f],digits),col=((sign(e.loadings[v,f])<0) +1),lty=((sign(e.loadings[v,f])<0)+1),adj=f %% adj +1)}
    			}
             }
   
   }
}  




#draw a heterarchy diagram 
"het.diagram" <- function(r,levels,cut=.3,digits=2,both=TRUE,main="Heterarchy diagram",l.cex ,gap.size,...) {
col = c("black","red")
nlevels <- length(levels)
if(missing(gap.size)) gap.size <- .4/nlevels
if(missing(l.cex)) l.cex <- 1
nvar <- max(unlist(lapply(levels,length)))
lowest <- rownames(r)[levels[[1]]]

xlim <- c(-.25,(nlevels-.75))
ylim <- c(.25,nvar)
plot(0,type="n",frame.plot=FALSE,axes=FALSE,xlim=xlim,ylim=ylim,ylab="",xlab="",main=main,...)
x <- 0

#first draw the left most layer
if(!is.null(names(levels))) {text(x,nvar,names(levels)[x+1]) }
nlevel.i <- length(levels[[1]])
lower <- list(nlevel.i)
spacing <- (nvar)/(nlevel.i+1)
for (y in 1:nlevel.i) {
lower[[y]] <- dia.rect(x,y*spacing,lowest[y],xlim=xlim,ylim=ylim,...)
}
names.i <- lowest
#now repeat for each higher layer
for(i in 2:(nlevels)){
nlevel.i <- length(levels[[i]])
level.i <- list(nlevel.i)
names.next <- rownames(r)[levels[[i]]]
x <- i-1
if(!is.null(names(levels))) {text(x,nvar,names(levels)[i]) }
spacing <- (nvar)/(nlevel.i+1)
for(y in 1:nlevel.i) {
	level.i[[y]] <- dia.rect(x,y*spacing,rownames(r)[levels[[i]][y]],xlim=xlim,ylim=ylim,...)	
    	for(j in 1:length(levels[[i-1]])) {
       		if(abs(r[names.i[j],names.next[y]]) > cut) { 
       		 dia.arrow(from=level.i[[y]]$left, to=lower[[j]]$right, labels=round(r[names.i[j],names.next[y]],digits),both=both,adj= y %%3 + 1,
       		col=(sign(r[names.i[j],names.next[y]] < 0) +1),lty=(sign(r[names.i[j],names.next[y]] < 0)+1),l.cex = l.cex,gap.size=gap.size,...)
   				} }
   } #end of drawing the factors for this  level
     	names.i <- names.next 
     	lower <-level.i			 
} #end of levels (i) loop
}
                