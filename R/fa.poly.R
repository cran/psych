 #polychoric factor analysis with confidence intervals
 "fa.poly" <- 
function(x,nfactors=1,n.obs = NA,n.iter=1,rotate="oblimin",SMC=TRUE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05,scores="regression",oblique.scores=TRUE,weight=NULL,global =TRUE, ...) {
 cl <- match.call()
 ncat <- 8
n.obs <- dim(x)[1]
	tx <- table(as.matrix(x))
	
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x)
	    typ = "tet"} else {
	    
	    tab <- apply(x,2,function(x) table(x))
if(is.list(tab)) {len <- lapply(tab,function(x) length(x))} else {len <- dim(tab)[1] }
nvar <- ncol(x)
dvars <- subset(1:nvar,len==2)   #find the dichotomous variables
pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)
if(length(pvars)==ncol(x)) {tet <- polychoric(x,weight=weight,global=global)
	    typ = "poly"} else {tet <- mixed.cor(x,weight=weight,global=global)
	    typ="mixed" }}
 r <- tet$rho
 #call fa with the polychoric/tetrachoric matrix
 #fa will not return scores, we still need to find them
 
 f <- fa(r,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC = SMC,missing=FALSE,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,scores=scores,oblique.scores=oblique.scores,...) #call fa with the appropriate parameters
 f$Call <- cl
 fl <- f$loadings  #this is the original
 nvar <- dim(fl)[1]

 if(n.iter > 1) {
 e.values <- list(pc =vector("list",n.iter),fa =vector("list",n.iter))
 replicates <- vector("list",n.iter)
 rep.rots <- vector("list",n.iter)
 for (trials in 1:n.iter) {
 xs <- x[sample(n.obs,n.obs,replace=TRUE),]
  if(typ!= "tet") {tets <- mixed.cor(xs,weight=weight,global=global)} else {tets <- tetrachoric(xs,weight=weight)}
  r <- tets$rho
  values.samp <- eigen(tets$rho)$values
   					e.values[["pc"]][[trials]] <- values.samp
   					
   				
  fs <- fa(r,nfactors=nfactors,rotate=rotate,SMC = SMC,missing=FALSE,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,...) #call fa with the appropriate parameters
     e.values[["fa"]][[trials]] <- fs$values
 if(nfactors > 1) {t.rot <- target.rot(fs$loadings,fl)
                  replicates[[trials]] <- t.rot$loadings
                   if(!is.null(fs$Phi)) {  phis <- fs$Phi  # should we rotate the simulated factor  correlations? 
                    rep.rots[[trials]] <- phis[lower.tri(phis)]}}  else {
                    replicates[[trials]] <- fs$loadings}
 }
 
replicates <- matrix(unlist(replicates),ncol=nfactors*nvar,byrow=TRUE)

if(!is.null( f$Phi) ) {rep.rots <- matrix(unlist(rep.rots),ncol=nfactors*(nfactors-1)/2,byrow=TRUE) 
       z.rot <- fisherz(rep.rots)
       means.rot <- colMeans(z.rot,na.rm=TRUE)
      sds.rot <- apply(z.rot,2,sd, na.rm=TRUE)
      sds.rot <- fisherz2r(sds.rot)
      ci.rot.lower <- means.rot + qnorm(p/2) * sds.rot
      ci.rot.upper <- means.rot + qnorm(1-p/2) * sds.rot
       means.rot <- fisherz2r(means.rot)
      ci.rot.lower <- fisherz2r(ci.rot.lower)
      ci.rot.upper <- fisherz2r(ci.rot.upper)
      
      ci.rot <- data.frame(lower=ci.rot.lower,upper=ci.rot.upper)
} else  {rep.rots <- NULL
         means.rot <- NULL
         sds.rot <- NULL
         z.rot <- NULL
         ci.rot <- NULL }



z.replicates <- fisherz(replicates)  #convert to z scores

means <- matrix(colMeans(z.replicates,na.rm=TRUE),ncol=nfactors)
sds <-  matrix(apply(z.replicates,2,sd,na.rm=TRUE),ncol=nfactors)

ci.lower <-  means + qnorm(p/2) * sds 
ci.upper <- means + qnorm(1-p/2) * sds
means <- fisherz2r(means)
sds <- fisherz2r(sds)
ci.lower <- fisherz2r(ci.lower)
ci.upper <- fisherz2r(ci.upper)
#ci.low.e <- apply(replicates,2, quantile, p/2)
#ci.up.e <- apply(replicates,2,quantile, (1-p/2))
#ci <- data.frame(lower = ci.lower, upper=ci.upper, low.e=ci.low.e, up.e=ci.up.e)
ci <- data.frame(lower = ci.lower,upper=ci.upper)
class(means) <- "loadings"
#class(sds) <- "loadings"

colnames(means) <- colnames(sds) <- colnames(fl)
rownames(means) <- rownames(sds) <- rownames(fl)

ei.pc <-describe(matrix(unlist(e.values$pc),ncol=nvar,byrow=TRUE))  #eigen values of pcs
ei.fa <- describe(matrix(unlist(e.values$fa),ncol=nvar,byrow=TRUE)) #eigen values of fa
e.stats <- list(ob.fa=f$values,ob.pc=f$e.values,pc=ei.pc,fa=ei.fa)
 

results <- list(fa = f,rho=tet$rho,tau=tet$tau,n.obs=n.obs,means = means,sds = sds,ci = ci, means.rot=means.rot,sds.rot=sds.rot,ci.rot=ci.rot,Call= cl,replicates=replicates,rep.rots=rep.rots,e.values=e.values,e.stats=e.stats)

class(results) <- c("psych","fa.ci")
} else {results <- list(fa = f,rho=r,tau=tet$tau,n.obs=n.obs,Call=cl)
       if(oblique.scores) {results$scores <- factor.scores(x,f=f$loadings,Phi=f$Phi,method=scores,rho=r) } else {results$scores <- factor.scores(x,f=f$Structure,method=scores,rho=r)}
       class(results) <- c("psych","fa")
       }
       
return(results)
 }
 #written May 3 2011
 #revised Sept 13, 2013 to allow for mixed cor input
 #and to find factor scores if data are given
 #corrected Sept 20, 2013 to do the ci on the fisher zs and then convert back to r

 
 