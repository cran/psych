"fa.multi" <- 
function(r,nfactors=3,nfact2=1,n.obs = NA,n.iter=1,rotate="oblimin",scores="regression", residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",...) {
 cl <- match.call()
 
if(nfactors < 2) stop("number of lower level factors must be at least 2")
f1 <- fa(r=r,nfactors=nfactors,n.obs=n.obs,rotate=rotate,scores=scores,residuals=residuals,SMC = SMC,covar=covar,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor, ...=...) #call fa with the appropriate parameters
f2 <- fa(f1$Phi,nfactors=nfact2,rotate=rotate,fm=fm)

result <- list(f1=f1,f2=f2)
return(result)
}



#based upon the omegaDiagram function

"fa.multi.diagram" <- function(multi.results,sort=TRUE,labels=NULL,flabels=NULL,cut=.2,gcut=.2,simple=TRUE,errors=FALSE,
    digits=1,e.size=.1,rsize=.15,side=3,main=NULL,cex=NULL,color.lines=TRUE
    ,marg=c(.5,.5,1.5,.5),adj=2, ...) {
     if(color.lines) { colors <- c("black","red")} else {colors <- c("black","black") }
  Phi <- NULL  #the default case
  if(is.null(cex)) cex <- 1
   old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back
    
 if(sort) multi.results$f1 <- fa.sort(multi.results$f1)   
factors <- as.matrix(multi.results$f1$loadings)
          if(is.null(main)) {main <- "Hierarchical (multilevel) Structure" }
       
       gloading <- multi.results$f2$loading
        nvar <- num.var <- dim(factors)[1]   #how many variables?
   num.factors <- dim(factors)[2]
   num.2level <- dim(gloading)[2]
    


   
   e.size <- e.size * 10/ nvar   #this is an arbitrary setting that seems to work
#first some basic setup parameters 
  
   
   
    vars <- paste("V",1:num.var,sep="")  
   if (!is.null(labels)) {vars <- paste(labels)} else{vars <- rownames(factors) }
   if(!is.null(flabels)) {fact <- flabels} else { fact <- c(paste("F",1:num.factors,sep="")) }   # e.g.  "g"  "F'1" "F2" "F3"
 
   colnames(factors)[1:length(fact)] <- fact
   var.rect <- list()
   fact.rect <- list()
    max.len <- max(nchar(rownames(factors)))*rsize
   cex <-  min(cex,40/nvar)
   xleft <-  -max.len/2 #0#
   xright <- nvar + 1  # or hard code to 3?
  plot(0,type="n",xlim=c(xleft-max.len,xright+1),ylim=c(1,nvar+1),frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
   
            vloc <- xleft
            gloc <- xright
            grouploc <- (xright)/2
            start <- 0
            end <- num.factors
          
  for (v in 1:nvar) { 
 	 var.rect[[v]] <- dia.rect(vloc,nvar-v+1,rownames(factors)[v],xlim=c(0,nvar),ylim=c(0,nvar),cex=cex,...)
     }
   f.scale <- (nvar+ 1)/(num.factors+1)
   f.shift <- nvar/num.factors
   for (f in 1:num.factors) {
   		fact.rect[[f]] <- dia.ellipse(grouploc,(num.factors+1-f)*f.scale,colnames(factors)[f+start],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
     		for (v in 1:nvar)  {
    			if (abs(factors[v,f+start]) > cut) {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,col=colors[((sign(factors[v,f+start])<0) +1)],lty=((sign(factors[v,f+start])<0)+1),labels=round(factors[v,f+start],digits),adj=f %% adj +1)
                               }
                             }
                     }

 g.ellipse <- list()
 f.2scale <- (num.var+1)/(num.2level +1)
 for(g in 1:num.2level) {  
  g.ellipse[[g] ]<-  dia.ellipse(gloc,(num.2level+1-g)*f.2scale,colnames(gloading)[g],xlim=c(0,nvar),ylim=c(0,nvar),e.size=e.size,...)
 
   for (f in 1:num.factors) {
        			if (abs(gloading[f,g]) > cut)  {         dia.arrow(from=g.ellipse[[g]],to=fact.rect[[f]],col=colors[((sign(gloading[f,g])<0) +1)],lty=((sign(gloading[f,g])<0) +1),labels=round(gloading[f,g],digits),adj=f %% adj +1)                
                            }}
                  }
 
  if (errors) {for (v in 1:nvar) {
       dia.self(location=var.rect[[v]],scale=.5,side=side)  }
       }
}                  

