"fa.sapa" <- 
function(r,nfactors=1,n.obs = NA,n.iter=1,rotate="oblimin",scores="regression", residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",correct=.5,weight=NULL,frac=.1,...) {
 cl <- match.call()
  if(isCorrelation(r)) {if(is.na(n.obs) && (n.iter >1)) stop("You must specify the number of subjects if giving a correlation matrix and doing confidence intervals")
                               #  if(!require(MASS)) stop("You must have MASS installed to simulate data from a correlation matrix")
                                 }
  
 f <- fac(r=r,nfactors=nfactors,n.obs=n.obs,rotate=rotate,scores=scores,residuals=residuals,SMC = SMC,covar=covar,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor, correct=.5,weight=weight,...=...) #call fa with the appropriate parameters
 fl <- f$loadings  #this is the original
 #f1 <- fa.sort(f1)   #put them into echelon form  But this does not work with target.rot
# if(!require(parallel)) {message("Parallels is required to do confidence intervals")}

 nvar <- dim(fl)[1]
 
 if(n.iter > 1) {
 if(is.na(n.obs) ) {n.obs <- f$n.obs} 
 replicates <- list()
 rep.rots <- list()
 
 
 #using cor="tet" seems to lead to an error being thrown in factoring, which in turn hangs mclapply
if(cor!="tet") {replicateslist <- parallel::mclapply(1:n.iter,function(x) {
 #replicateslist <- lapply(1:n.iter,function(x) {
 if(isCorrelation(r)) {#create data sampled from multivariate normal with observed correlation
                                      mu <- rep(0, nvar)
                                      #X <- mvrnorm(n = n.obs, mu, Sigma = r, tol = 1e-06, empirical = FALSE)
                                      #the next 3 lines replaces mvrnorm (taken from mvrnorm, but without the checks)
                                      eX <- eigen(r)
                                      X <- matrix(rnorm(nvar * n.obs),n.obs)
                                      X <-  t(eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(X))
                            } else { X <- r[sample(n.obs,n.obs*frac,replace=FALSE),]}
  fs <- fac(X,nfactors=nfactors,rotate=rotate,scores="none",SMC = SMC,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor,correct=correct,...=...) #call fa with the appropriate parameters
  if(nfactors == 1) {  npairs <- pairwiseCount(X,diagonal=FALSE)
                 mean.npairs <- mean(npairs,na.rm=TRUE)
                 replicates <- list(loadings=fs$loadings,npairs=mean.npairs)
  } else  {
                    t.rot <- target.rot(fs$loadings,fl)
                     npairs <- pairwiseCount(X,diagonal=FALSE)
                    mean.npairs <- mean(npairs,na.rm=TRUE)
                   if(!is.null(fs$Phi)) {  phis <- fs$Phi  # should we rotate the simulated factor  correlations?
                   #we should report the target rotated phis, not the untarget rotated phis 
                     replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(t.rot$Phi)],npairs=mean.npairs)   #corrected 6/10/15
                    #replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(phis)])
                    }  else 
                   {replicates <- list(loadings=t.rot$loadings)}                
  }}) } else {replicateslist <- lapply(1:n.iter,function(x) {  #avoiding parallel for this case
 if(isCorrelation(r)) {#create data sampled from multivariate normal with observed correlation
                                      mu <- rep(0, nvar)
                                      #X <- mvrnorm(n = n.obs, mu, Sigma = r, tol = 1e-06, empirical = FALSE)
                                      #the next 3 lines replaces mvrnorm (taken from mvrnorm, but without the checks)
                                      eX <- eigen(r)
                                      X <- matrix(rnorm(nvar * n.obs),n.obs)
                                      X <-  t(eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(X))
                            } else {X <- r[sample(n.obs,n.obs*frac,replace=FALSE),]}
  fs <- fac(X,nfactors=nfactors,rotate=rotate,scores="none",SMC = SMC,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor,correct=correct,...=...) #call fa with the appropriate parameters
  if(nfactors == 1) {replicates <- list(loadings=fs$loadings)} else  {
                    t.rot <- target.rot(fs$loadings,fl)
                
                 npairs <- pairwiseCount(X,diagonal=FALSE)
                 mean.npairs <- mean(npairs,na.rm=TRUE)
                   if(!is.null(fs$Phi)) {  phis <- fs$Phi  # should we rotate the simulated factor  correlations?
                   #we should report the target rotated phis, not the untarget rotated phis 
                     replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(t.rot$Phi)],npairs=mean.npairs)   #corrected 6/10/15
                    #replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(phis)])
                    }  else 
                   {replicates <- list(loadings=t.rot$loadings)}                
  }})}
  

replicates <- matrix(unlist(replicateslist),nrow=n.iter,byrow=TRUE)
mean.npairs <- mean(replicates[,NCOL(replicates)],na.rm=TRUE)
sds.pairs <- sd(replicates[,NCOL(replicates)],na.rm=TRUE)
replicates <- replicates[,-NCOL(replicates),drop=FALSE]


#drop weird replications (loadings > 1)
replicates[abs(replicates) > 1]  <- NA
means <- colMeans(replicates,na.rm=TRUE)
sds <- apply(replicates,2,sd,na.rm=TRUE)

if(length(means) > (nvar * nfactors ) ) {
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
   
   means <- matrix(means[1:(nvar*nfactors)],ncol=nfactors)
   sds <- matrix(sds[1:(nvar*nfactors)],ncol=nfactors)
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
f$cis <- list(means = means,sds = sds,ci = ci,p =2*ptci, means.rot=means.rot,sds.rot=sds.rot,ci.rot=ci.rot,p.rot = ptcirot,Call= cl,replicates=replicates,rep.rots=rep.rots,mean.pair=mean.npairs,sds.pairs=sds.pairs)
results <- f 
 results$Call <- cl
class(results) <- c("psych","fa.ci")
} else {results <- f
        results$Call <- cl
       class(results) <- c("psych","fa")
       }
return(results)

 }
 #written May 1 2011
 #modified May 8, 2014 to make cis an object in f to make sorting easier
 