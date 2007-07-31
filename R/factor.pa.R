"factor.pa" <- 
function(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NULL, min.err = .001,digits=2,max.iter=50) {
    n <- dim(r)[2]
    if (n!=dim(r)[1]) {  n.obs <- dim(r)[1]
    		r <- cor(r,use="pairwise") # if given a rectangular matrix, the find the correlations first
           }
     if(!is.matrix(r)) {r <- as.matrix(r)}
    if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)}

    orig <- diag(r)
    r.mat <- r
    comm <- sum(diag(r.mat))
    err <- comm
     i <- 1
    comm.list <- list()
    while(err > min.err)    #iteratively replace the diagonal with our revised communality estimate
      {
        eigens <- eigen(r.mat)
        #loadings <- eigen.loadings(eigens)[,1:nfactors]
         if(nfactors >1 ) {loadings <- eigens$vectors[,1:nfactors] %*% diag(sqrt(eigens$values[1:nfactors])) } else {loadings <- eigens$vectors[,1] * sqrt(eigens$values[1] ) }
         model <- loadings %*% t(loadings)
         
         new <- diag(model)
         
         comm1 <- sum(new)
         diag(r.mat) <- new
         err <- abs(comm-comm1)
         comm <- comm1
         comm.list[[i]] <- comm1
         i <- i + 1
         if(i > max.iter) {warning("maximum iteration exceeded")
                err <-0 }
       }

       #make each vector signed so that the maximum loading is positive
    if (nfactors >1) {
    		maxabs <- apply(apply(loadings,2,abs),2,which.max)
   			sign.max <- vector(mode="numeric",length=nfactors)
    		for (i in 1: nfactors) {sign.max[i] <- sign(loadings[maxabs[i],i])}
    		loadings <- loadings %*% diag(sign.max)
   		
    	} else {
    		mini <- min(loadings)
   			maxi <- max(loadings)
    		if (abs(mini) > maxi) {loadings <- -loadings }
    		loadings <- as.matrix(loadings)
    	} #sign of largest loading is positive
    	colnames(loadings) <- paste("PA",1:nfactors,sep='')
    rownames(loadings) <- rownames(r)
    
    model <- loadings %*% t(loadings)  
    result$communality <- round(diag(model),digits)
    
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
    if(nfactors<1) nfactors <- n
   
    residual<- r - model
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
    if (residuals) {result$residual <- round(residual,digits)}
  
    r2.off <- r
    diag(r2.off) <- 0
    r2.off <- sum(r2.off^2)
    diag(residual) <- 0
    rstar.off <- sum(residual^2)
    result$fit <- round(1-rstar2/r2,digits)
    result$fit.off <- round(1-rstar.off/r2.off,digits)
    result$values <- round(eigens$values,digits)
  
    diag(model) <- 1   
    model.inv <- solve(model)
    m.inv.r <- model.inv %*% r
    result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n 
    if (!is.null(n.obs)) {result$STATISTIC <-  result$objective * (n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3
    	if (result$dof > 0) {result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else result$PVAL <- NA}
    result$loadings <- round(loadings,digits)
    if (rotate =="oblimin") {result$phi <- phi}
    result$communality.iterations <- round(unlist(comm.list),digits)
    return(result) }