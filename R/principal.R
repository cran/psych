"principal" <-
function(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NA, scores=FALSE,missing=FALSE,impute="median", digits=2) {
   n <- dim(r)[2]
  
   if (n!=dim(r)[1]) {
      				n.obs <- dim(r)[1] 
      				if(scores) {x.matrix <- r
    
                  if(missing) {
        			miss <- which(is.na(x.matrix),arr.ind=TRUE)
       				 if(impute=="mean") {
       				item.means <- colMeans(x.matrix,na.rm=TRUE)   #replace missing values with means
      				x.matrix[miss]<- item.means[miss[,2]]} else {
       				item.med   <- apply(x.matrix,2,median,na.rm=TRUE) #replace missing with medians
       				 x.matrix[miss]<- item.med[miss[,2]]}
       				 }}
      r <- cor(r,use="pairwise")  # if given a rectangular matrix, then find the correlations first
             } else {
     				if(!is.matrix(r)) {  r <- as.matrix(r)}
     				 sds <- sqrt(diag(r))    #convert covariance matrices to correlation matrices
                     r <- r/(sds %o% sds)  }  #added June 9, 2008
    if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0,fit.off=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0,fit.off=0)}
    eigens <- eigen(r)    #call the eigen value decomposition routine
    result$values <- round(eigens$values,digits)
    loadings <- eigens$vectors %*% sqrt(diag(eigens$values))
   if(nfactors >0) {loadings <- loadings[,1:nfactors]} else {nfactors <- n}

   	
   
    	 
     #added January 19, 2009 to flip based upon colSums of loadings
    if (nfactors >1) {sign.tot <- vector(mode="numeric",length=nfactors)
                 sign.tot <- sign(colSums(loadings))
                 loadings <- loadings %*% diag(sign.tot)
     } else { if (sum(loadings) <0) {loadings <- -as.matrix(loadings)} else {loadings <- as.matrix(loadings)}
             colnames(loadings) <- "PC1" }
     
    	
    colnames(loadings) <- paste("PC",1:nfactors,sep='')
    rownames(loadings) <- rownames(r)

     Phi <- NULL
     if(rotate != "none") {if (nfactors > 1) {
   			if (rotate=="varimax") { 
   			loadings <- varimax(loadings)$loadings 
   			 Phi <- NULL} else { 
     			if ((rotate=="promax") | (rotate=="Promax"))  {pro <- Promax(loadings)
     			                loadings <- pro$loadings
     			                 Phi<-pro$Phi 
     			                  } else {
     			if (rotate =="oblimin") {
     				if (!require(GPArotation)) {warning("I am sorry, to do oblimin rotations requires the GPArotation package to be installed")
     				Phi <- NULL} else { ob  <- oblimin(loadings)
     				loadings <- ob$loadings
     				 Phi <- ob$Phi}
     		                             }
     	               }}
     			
     }}
   #just in case the rotation changes the order of the components, sort them by size of eigen value
   if(nfactors >1) {
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated,decreasing=TRUE)
    loadings <- loadings[,ev.order]}
    if(!is.null(Phi)) {Phi <- Phi[ev.order,ev.order] } #January 20, 2009 but, then, we also need to change the order of the rotation matrix!
   
     class(loadings) <- "loadings"
     
    model <- loadings %*% t(loadings)  
    result$communality <- diag(model)
    result$uniquenesses <- diag(r- model)
    residual<- r - model
    if (residuals) {result$residual <-residual}
    diag(model) <- 1   
    model.inv <- solve(model)
    m.inv.r <- model.inv %*% r
    if(is.na(n.obs)) {result$n.obs=NA
    			        result$PVAL=NA} else {result$n.obs=n.obs}
    result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
  
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n 
    result$criteria <- c("objective"=result$objective,NA,NA)
    if (!is.na(n.obs)) {result$STATISTIC <-  result$objective * ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3)
            if (result$STATISTIC <0) {result$STATISTIC <- 0}  
   			if (result$dof > 0) {result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else result$PVAL <- NA}
     
      if(!is.null(Phi)) {result$Phi <- phi}
  
    
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
   diag(residual) <- 0
   rstar.off <- sum(residual*residual)
   r2.off <- r2 - n
   
    result$fit <- round(1-rstar2/r2,digits)
    result$fit.off <- round(1-rstar.off/r2.off,digits)
    result$loadings <- round(loadings,digits)
    if (!is.null(Phi)) {result$Phi <- Phi}
    result$factors <- nfactors 
    class(result) <- c("psych", "principal") 
    result$fn <- "principal"
    if(scores) {result$scores <- factor.scores(scale(x.matrix),loadings) }
    return(result)
   }
  
  # modified August 10, 2007
  # modified Feb 2, 2008 to calculate scores and become a factanal class
  #Modified June 8,2008 to get chi square values to work or default to statistic if n.obs==NULL.
  #modified Jan  2, 2009 to  report the correlations between oblique factors
  #modified December 30 to let n.obs ==NA rather than NULL to be compatible with factanal
  #modified Jan 14, 2009 to change phi to Phi  to avoid confusion
  #modified Jan 19, 2009 to allow for promax rotations