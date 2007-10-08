"principal" <-
function(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NULL, digits=2) {
   n <- dim(r)[2]
  
   if (n!=dim(r)[1]) {
      n.obs <- dim(r)[1] 
      r <- cor(r,use="pairwise") } # if given a rectangular matrix, then find the correlations first
    
    if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0,fit.off=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0,fit.off=0)}
    eigens <- eigen(r)    #call the eigen value decomposition routine
    result$values <- round(eigens$values,digits)
    loadings <- eigens$vectors %*% sqrt(diag(eigens$values))
   if(nfactors >0) {loadings <- loadings[,1:nfactors]} else {nfactors <- n}

   	
   	if(nfactors >1) {
    	maxabs <- apply(apply(loadings,2,abs),2,which.max)
    	sign.max <- vector(mode="numeric",length=nfactors)
    	for (i in 1: nfactors) {sign.max[i] <- sign(loadings[maxabs[i],i])}
  
  		loadings <- loadings %*% diag(sign.max)  #sign of largest loading is positive 
      	} else {if(nfactors >0) {mini <- min(loadings)
   			maxi <- max(loadings)
    		if (abs(mini) > maxi) {loadings <- -loadings }
    		loadings <- as.matrix(loadings)
    	 } }
    	
    colnames(loadings) <- paste("PC",1:nfactors,sep='')
    rownames(loadings) <- rownames(r)


     if(rotate != "none") {if (nfactors > 1) {
    	if (rotate=="varimax") { 
   			loadings <- varimax(loadings)$loadings } else { 
     			if (rotate=="promax") {loadings <- promax(loadings)$loadings } else {
     			if (rotate =="oblimin") {
     				if (!require(GPArotation)) {warning("I am sorry, to do oblimin rotations requires the GPArotation package to be installed")
     				phi <- NULL} else { ob  <- oblimin(loadings)
     				 loadings <- ob$loadings
     				 phi <- ob$Phi}
     		                             }
     	               }}
     }}
   
    class(loadings) <- "loadings"
    model <- loadings %*% t(loadings)  
    result$communality <- diag(model)
    residual<- r - model
    if (residuals) {result$residual <-residual}
    diag(model) <- 1   
    model.inv <- solve(model)
    m.inv.r <- model.inv %*% r
    result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n 
    if (!is.null(n.obs)) {result$STATISTIC <-  result$objective * (n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3
   if (result$dof >0) { result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else result$PVAL  < - NA}
    
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
   diag(residual) <- 0
   rstar.off <- sum(residual*residual)
   r2.off <- r2 - n
   
    result$fit <- round(1-rstar2/r2,digits)
    result$fit.off <- round(1-rstar.off/r2.off,digits)
    result$loadings <- round(loadings,digits)
    if (rotate =="oblimin") {result$phi <- phi}
    return(result)
   }
  
  #last modified August 10, 2007 
