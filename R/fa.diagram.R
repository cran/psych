#factor analysis and sem diagrams
#based upon fa.graph with some ideas taken from the diagram and shape packages of Karline Soetaert 
#version of September 30, 2009
#developed to replace Rgraphviz which is too much of a pain to install
#Rgraphviz uses a NEL  (Node, Edge, Label) representation while diagram uses a complete linking matrix
#thus, I am trying to combine these two approaches
#some small revisions, March 2015 to allow cex to be passed to dia.ellipse and dia.rect

"fa.diagram" <-
  function(fa.results,Phi=NULL,fe.results=NULL,sort=TRUE,labels=NULL,cut=.3,simple=TRUE,errors=FALSE,g=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main,cex=NULL,marg=c(.5,.5,1,.5),adj=1, ...) {
    old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back
   col <- c("black","red")
   if(missing(main)) if(is.null(fe.results)) {main <- "Factor Analysis" } else {main <- "Factor analysis and extension"}
   if(!is.matrix(fa.results) && !is.null(fa.results$fa) && is.list(fa.results$fa)) fa.results <- fa.results$fa
   if(is.null(cex)) cex <- 1
  #Phi <- NULL  #the default case
 if(sort) {fa.results <- fa.sort(fa.results) 
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
 	 var.rect[[v]] <- dia.rect(left,top -v - max(0,n.evar-nvar)/2  ,rownames(factors)[v],xlim=limx,ylim=limy,cex=cex,...)
     }
   f.scale <- (top)/(num.factors+1)
   f.shift <- max(nvar,n.evar)/num.factors
   if(g) {fact.rect[[1]] <- dia.ellipse(-max.len/2,top/2,colnames(factors)[1],xlim=limx,ylim=limy,e.size=e.size,cex=cex,...)
          	for (v in 1:nvar)  {if(simple && (abs(factors[v,1]) == max(abs(factors[v,])) )  && (abs(factors[v,1]) > cut) | (!simple && (abs(factors[v,1]) > cut))) { 
    			dia.arrow(from=fact.rect[[1]],to=var.rect[[v]]$left,labels =round(factors[v,1],digits),col=((sign(factors[v,1])<0) +1),lty=((sign(factors[v,1])<0)+1))
    	 }}}
   for (f in gf:num.factors) {
   		fact.rect[[f]] <- dia.ellipse(left+middle,(num.factors+gf-f)*f.scale,colnames(factors)[f],xlim=limx,ylim=limy,e.size=e.size,cex=cex,...)
     		for (v in 1:nvar)  {
     		
    			if(simple && (abs(factors[v,f]) == max(abs(factors[v,])) )  && (abs(factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) { 
    			dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =round(factors[v,f],digits),col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0)+1),adj=f %% adj +1,cex=cex)
    			 }
   }
   }
   
   if(!is.null(Phi) && (ncol(Phi) >1)) { for (i in 2:num.factors) {
     for (j in 1:(i-1)) {
     if(abs(Phi[i,j]) > cut) {
       # dia.curve(from=c(x.max-2+ e.size*nvar,(num.factors+1-i)*f.scale),to=c(x.max -2+ e.size*nvar,(num.factors+1-j)*f.scale),labels=round(Phi[i,j],digits),scale=(i-j),...)}
		dia.curve(from=fact.rect[[j]]$right,to=fact.rect[[i]]$right,labels=round(Phi[i,j],digits),scale=(i-j),cex=cex,...)} }
  															 }
 
						}
	
  if (errors) {for (v in 1:nvar) {
       dia.self(location=var.rect[[v]],scale=.5,side=side)  }
       }
       
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
                