 #Pearson or polychoric correlations with confidence intervals

 "cor.ci" <- 
function(x, keys = NULL, n.iter = 100,  p = 0.05, poly = FALSE, method = "pearson") {
 cl <- match.call()
 n.obs <- dim(x)[1]
 
 
if(poly) {
 ncat <- 8
nvar <- dim(x)[2]
	tx <- table(as.matrix(x))
	
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x)
	    typ = "tet"} else {
	    
	    tab <- apply(x,2,function(x) table(x))
       if(is.list(tab)) {len <- lapply(tab,function(x) length(x))} else {len <- dim(tab)[1] }

dvars <- subset(1:nvar,len==2)   #find the dichotomous variables
pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)
if(length(pvars)==ncol(x)) {tet <- polychoric(x)
	    typ = "poly"} else {tet <- mixed.cor(x)
	    typ="mixed" }
	    
	  } 
   
    rho <- tet$rho 
 } else { rho <- cor(x,use="pairwise",method=method)
 
         } 
         if(!is.null(keys)) {bad <- FALSE
         if(any(is.na(rho))) {warning("Some of the item correlations are NA and thus finding scales that include those items will not work.\n We will try to do it for those  scales which do not include those items.
         \n Proceed with caution. ") 
         bad <- TRUE
         rho <- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(rho*x))*x,na.rm=TRUE))  #matrix multiplication without matrices!
         } else {
          rho <- t(keys) %*% rho %*% keys} }  #find the correlation between the scales

 rho <- cov2cor(rho)     #scale to correlations  
  nvar <- dim(rho)[2]
 
 if(n.iter > 1) {
 
 replicates <- list()
 rep.rots <- list()
 
if(!require(parallel)) {warning("parallel package needed for mclapply")}
 replicates <- mclapply(1:n.iter,function(XX) {
 xs <- x[sample(n.obs,n.obs,replace=TRUE),]
 {if(poly) {
  if(typ!= "tet") {tets <- mixed.cor(xs)} else {tets <- tetrachoric(xs)}

  R <- tets$rho} else {R <- cor(xs,use="pairwise",method=method)}
 
 if(!is.null(keys)) { if (bad) {covariances <- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(R*x))*x,na.rm=TRUE))  #matrix multiplication without matrices!
  } else {
 covariances <- t(keys) %*% R %*% keys}
               r <- cov2cor(covariances) } else {r <- R}					
             
 rep.rots <- r[lower.tri(r)]
 }
 } 
 )
 
 }

   replicates <- matrix(unlist(replicates),ncol=nvar*(nvar-1)/2,byrow=TRUE) 
       z.rot <- fisherz(replicates)
     means.rot <- colMeans(z.rot,na.rm=TRUE)
      sds.rot <- apply(z.rot,2,sd, na.rm=TRUE)
      sds.rot <- fisherz2r(sds.rot)
      ci.rot.lower <- means.rot + qnorm(p/2) * sds.rot
      ci.rot.upper <- means.rot + qnorm(1-p/2) * sds.rot
      means.rot <- fisherz2r(means.rot)
      ci.rot.lower <- fisherz2r(ci.rot.lower)
      ci.rot.upper <- fisherz2r(ci.rot.upper)
      low.e <- apply(replicates,2,quantile, p/2,na.rm=TRUE)
      up.e  <- apply(replicates, 2, quantile, 1-p/2,na.rm=TRUE)
      ci.rot <- data.frame(lower=ci.rot.lower,low.e=low.e,upper=ci.rot.upper,up.e=up.e)
      cnR <- abbreviate(colnames(rho),minlength=5) 
      k <- 1
     for(i in 1:(nvar-1)) {for (j in (i+1):nvar) {
      rownames(ci.rot)[k] <- paste(cnR[i],cnR[j],sep="-")
      k<- k +1 }}
      
     


results <- list(rho=rho, means=means.rot,sds=sds.rot,ci=ci.rot,Call= cl,replicates=rep.rots)

class(results) <- c("psych","cor.ci")
    
return(results)

 }
 #written Sept 20, 2013 
 #adapted from fa.poly
 
 

 