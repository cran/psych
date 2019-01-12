

"fa.extension" <-
  function(Roe,fo,correct=TRUE) {
 cl <- match.call()
 omega <-FALSE

 if(!is.null(class(fo)[2])) {if(class(fo)[2]=="fa") {
          if(!is.null(fo$Phi)) {Phi <- fo$Phi} else {Phi <- NULL}
          
       fl <- fo$loadings 
       fs <- fo$Structure
       
     } else {if (class(fo)[2] == "omega") {
         omega <- TRUE
         w <- fo$stats$weights
         fl <- fo$schmid$sl
         Phi <- NULL
         fl <- fl[,1:(dim(fl)[2]-3)]
         nfactors <- dim(fl)[2]
         fe <- t(t(w) %*% Roe)
         foblique <- fo$schmid$oblique
         feoblique <- t( Roe) %*% foblique %*% (solve(t(foblique)%*% (foblique)))
         feoblique <- feoblique %*% solve(fo$schmid$phi)  
       } 
    }
    }
 
 if(!omega) fe <- t( Roe) %*% fl %*% (solve(t(fl)%*% (fl))) 
  if(!is.null(Phi)) fe <- fe %*% solve(Phi)
 
 if(!correct) {#the Gorsuch case
     d <-diag(t(fl) %*% fo$weight)
     fe <- (fe * d)
 }
 colnames(fe) <- colnames(fl)
 rownames(fe) <- colnames(Roe)
if(!is.null(Phi)) {resid <- Roe - fl %*% Phi %*% t(fe)} else {resid <- Roe - fl  %*% t(fe)}   #fixed to actually give residual  (1/30/18)
 result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)
 if(!omega)  {result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)} else {result <- list(loadings = fe,oblique= feoblique,Phi=Phi,resid=resid,Call=cl)}
 class(result) <- c("psych","extension")
 return(result)
}
#written April 5, 2011
#revised August 15, 2011 to avoid using the weights matrix except in the omega case

#created December 8, 2012 to allow for extension and goodness of fits of total model
#modified 31/5/14 to allow for omega extension as well 
#modified 04-09/16 to pass the Structure matrix as well
"fa.extend" <- 
function(r,nfactors=1,ov=NULL,ev=NULL,n.obs = NA, np.obs=NULL,correct=TRUE,rotate="oblimin",SMC=TRUE,warnings=TRUE, fm="minres",alpha=.1, omega=FALSE, ...) {
 cl <- match.call()
  nv <- c(ov,ev)
 if(nrow(r) > ncol(r)){  #the case of a data matrix
 if(omega) {fo <- omega(r[,ov],nfactors=nfactors,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,...)} else {
       fo <- fa(r[,ov],nfactors=nfactors,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,...)}
         n.obs <- nrow(r)
         np.obs.r <- pairwiseCount(r)[nv,nv]
         np.obs <- np.obs.r[ov,ov]
        r <- cor(r,use='pairwise')  
    } else {  #the case of a correlation matrix         
       R <- r[ov,ov]
       np.obs.r <- np.obs
      if(omega) {fo <- omega(R,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,np.obs=np.obs[ov,ov],...)} else { 
      fo <- fa(R,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,np.obs=np.obs[ov,ov],...)}
     }
Roe <- r[ov,ev,drop=FALSE]
fe <- fa.extension(Roe,fo,correct=correct)
if(omega) fo$loadings <- fo$schmid$sl[,1:(ncol(fo$schmid$sl)-3)]
       
foe <- rbind(fo$loadings,fe$loadings)

if(omega) oblique <- rbind(fo$schmid$oblique,fe$oblique)

if(is.na(n.obs) && !is.null(np.obs)) n.obs <- max(as.vector(np.obs))
result <- factor.stats(r[nv,nv],foe,fo$Phi,n.obs,np.obs.r,alpha=alpha)
if(omega) result$schmid$sl <- foe
    result$rotation <- rotate
    result$loadings <- foe
    if(nfactors > 1) {if(is.null(fo$Phi)) {h2 <- rowSums(foe^2)} else {h2 <- diag(foe %*% fo$Phi %*% t(foe)) }} else {h2 <-foe^2}
    result$communality <- h2
    result$fm <- fm  #remember what kind of analysis we did
    
    result$fo=fo
    if(omega) {result$schmid$sl <- foe
              result$schmid$gloading <- fo$schmid$gloading
              result$schmid$oblique <- oblique
            }
   if(is.null(fo$Phi)) {result$Structure <- foe } else { result$Structure <- foe %*% fo$Phi}
    result$fe=fe
    result$resid=fe$resid
    result$Phi=fo$Phi
    result$fn="fa"
    result$Call=cl
class(result) <- c("psych","extend")
return(result)
}


#adapted from fa.diagram but treats the extension variables as y variables
#draw the standard fa.diagram for the original variables and then regressions to the fe variables
#basically for the case of extension to criterion variables with lower path strengths 
#offers a bit more control in the e.cut and e.simple options

"extension.diagram" <-
  function(fa.results,Phi=NULL,fe.results=NULL,sort=TRUE,labels=NULL,cut=.3,e.cut=.1,simple=TRUE,e.simple=FALSE,errors=FALSE,g=FALSE,
    digits=1,e.size=.05,rsize=.15,side=2,main,cex=NULL,marg=c(.5,.5,1,.5),adj=1,ic=FALSE, ...) {
   # if(length(class(fa.results)) > 1) {if(class(fa.results)[2] == 'principal') {pc <- TRUE} else {pc <- FALSE}} else { pc <- FALSE}
   # if(ic) pc <- TRUE
   pc <- FALSE
    old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back
   col <- c("black","red")
   if(missing(main)) {main <- "Factor analysis and extension"}
 #  if(!is.matrix(fa.results) && !is.null(fa.results$fa) && is.list(fa.results$fa)) fa.results <- fa.results$fa
   if(is.null(cex)) cex <- 1
  #Phi <- NULL  #the default case
 if(sort) {
          
         fe.results <- fa.sort(fa.results$fo)} 
 if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {factors <- as.matrix(fe.results$loadings)
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
   for (f in gf:num.factors) {  #body  34
   		if (pc) {fact.rect[[f]] <- dia.rect(left+middle,(num.factors+gf-f)*f.scale,colnames(factors)[f],xlim=limx,ylim=limy,cex=cex,...) 
   		} else {fact.rect[[f]] <- dia.ellipse(left+middle,(num.factors+gf-f)*f.scale,colnames(factors)[f],xlim=limx,ylim=limy,e.size=e.size,cex=cex,...)}
     		for (v in 1:nvar)  {
     		
    			if(simple && (abs(factors[v,f]) == max(abs(factors[v,])) )  && (abs(factors[v,f]) > cut) | (!simple && (abs(factors[v,f]) > cut))) { 
    		if(pc) {dia.arrow(to=fact.rect[[f]],from=var.rect[[v]]$right,labels =round(factors[v,f],digits),col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0)+1),adj=f %% adj ,cex=cex) 
    		} else {dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$right,labels =round(factors[v,f],digits),col=((sign(factors[v,f])<0) +1),lty=((sign(factors[v,f])<0)+1),adj=f %% adj +1,cex=cex)}
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
     e.loadings <- fa.results$fe$loadings
    n.evar <- NROW(e.loadings)
    cut <- e.cut    #draw all extensions
    simple <- e.simple 
       
     for (v in 1:n.evar) { 
 	 var.rect[[v]] <- dia.rect(x.max,v* nvar/(n.evar+1),rownames(e.loadings)[v],xlim=limx,ylim=limy,cex=cex,...)
 	 for(f in 1:num.factors) {
 	 if(simple && (abs(e.loadings[v,f]) == max(abs(e.loadings[v,])) )  && (abs(e.loadings[v,f]) > cut) | (!simple && (abs(e.loadings[v,f]) > cut))) { 
    			dia.arrow(from=fact.rect[[f]],to=var.rect[[v]]$left,labels =round(e.loadings[v,f],digits),col=((sign(e.loadings[v,f])<0) +1),lty=((sign(e.loadings[v,f])<0)+1),adj=f %% adj +1)}
    			}
             }
   
   }
}  
 
 
