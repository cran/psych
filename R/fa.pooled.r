fa.pooled <- function(datasets,nfactors=1,rotate="oblimin",scores="regression", residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",correct=.5,weight=NULL,...) {

 cl <- match.call()
 replicates <- list()
 replicateslist <- list()
 rep.rots <- list()
 n.iter <- length(datasets)

#the first fa becomes the target for the remaining ones 
 X <- datasets[[1]]
 f <- fac(X,nfactors=nfactors,rotate=rotate,scores="none",SMC = SMC,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor,correct=correct,...=...) 
 
 fl <- f$loadings
nvar <- ncol(X)


 #now do the replicated
 for (iter in (1:n.iter)) {
 X <- datasets[[iter]]
fs <-  fac(X,nfactors=nfactors,rotate=rotate,scores="none",SMC = SMC,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor,correct=correct,...=...) #call fa with the appropriate parameters

 if(nfactors == 1) {replicateslist[[iter]] <- list(loadings=fs$loadings)} else  {
                    t.rot <- target.rot(fs$loadings,fl)
                
                   if(!is.null(fs$Phi)) {  phis <- fs$Phi  # should we rotate the simulated factor  correlations?
                   #we should report the target rotated phis, not the untarget rotated phis 
                     replicateslist[[iter]] <- list(loadings=t.rot$loadings,phis=phis[lower.tri(t.rot$Phi)])   #corrected 6/10/15
                    #replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(phis)])
                    }  else 
                   {replicateslist[[iter]] <- list(loadings=t.rot$loadings)}                
      t.rotstr <- target.rot(fs$Structure,fstr)
          Structure = t.rotstr$loadings
  }
  col_means <- colMeans(replicates, na.rm = TRUE)
  col_sds <- apply(replicates, 2, sd, na.rm = TRUE)
}

replicates <- matrix(unlist(replicateslist),nrow=n.iter,byrow=TRUE)

means <- colMeans(replicates,na.rm=TRUE)
sds <- apply(replicates,2,sd,na.rm=TRUE)

if(length(means) > (nvar * nfactors) ) {
   	means.rot <- means[(nvar*nfactors +1):length(means)]
   	sds.rot <-      sds[(nvar*nfactors +1):length(means)]  
	ci.rot.lower <- means.rot + qnorm(p/2) * sds.rot
  	ci.rot.upper <- means.rot + qnorm(1-p/2) * sds.rot  
   	ci.rot <- data.frame(lower=ci.rot.lower,upper=ci.rot.upper)    } else  {
        rep.rots <- NULL
        means.rot <- NULL
        sds.rot <- NULL
        z.rot <- NULL
        ci.rot <- NULL }
   
   tci <- abs(means)/sds
    ptci <- 1-pnorm(tci)
    if(!is.null(rep.rots)) {
   tcirot <- abs(means.rot)/sds.rot
   ptcirot <- 1- pnorm(tcirot)} else  {tcirot <- NULL
                                      ptcirot <- NULL}
ci.lower <-  means + qnorm(p/2) * sds
ci.upper <- means + qnorm(1-p/2) * sds

ci <- data.frame(lower = ci.lower,upper=ci.upper)
class(means) <- "loadings"

colnames(means) <- colnames(sds) <- colnames(fl)
rownames(means) <- rownames(sds) <- rownames(fl)

f$cis <- list(means = means,sds = sds,ci = ci,p =2*ptci, means.rot=means.rot,sds.rot=sds.rot,ci.rot=ci.rot,p.rot = ptcirot,Call= cl,replicates=replicates,rep.rots=rep.rots)
results <- f 
 results$Call <- cl
class(results) <- c("psych","fa.ci")

return(results)
}

