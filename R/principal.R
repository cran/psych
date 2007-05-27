"principal" <-
function(r,nfactors=0,residuals=FALSE,rotate=FALSE,digits=2) {
   n <- dim(r)[2]
  
   if (n!=dim(r)[1]) r <- cor(r,use="pairwise") # if given a rectangular matrix, then find the correlations first
    
    if (!residuals) { result <- list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)}
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


    if(rotate) {loadings <- varimax(loadings)$loadings }  
   
    residual<- factor.residuals(r,loadings)
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
    if (residuals) {result$residual <-residual}
    result$fit <- round(1-rstar2/r2,digits)
    result$loadings <- round(loadings,digits)
    return(result)
   }
  
  #last modified May 23, 2007 
