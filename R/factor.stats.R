"factor.stats" <- 
function(r,f,phi=NULL,n.obs=NA) {
cl <- match.call()
 n <- dim(r)[2]  #number of variables
 if(dim(r)[1] !=n ) {n.obs = dim(r)[1]
                    r <- cor(r,use="pairwise")
                     }
 if(is.data.frame(r)) r<- as.matrix(r)
 if((!is.matrix(f)) && (!is.data.frame(f)))  {f <- as.matrix(f$loadings)} else {f <- as.matrix(f)}
 nfactors <- dim(f)[2]  # number of factors
 if(is.null(phi)) {model <- f %*%  t(f)} else {model <- f %*% phi %*% t(f)}

 residual<- r - model
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
    result <- list(residual =residual)
  
    r2.off <- r
    diag(r2.off) <- 0
    r2.off <- sum(r2.off^2)
    diag(residual) <- 0
    rstar.off <- sum(residual^2)
    result$fit <-1-rstar2/r2
    result$fit.off <- 1-rstar.off/r2.off
    
    result$factors <- nfactors
  
    diag(model) <- 1   
     model.inv <- try(solve(model),silent=TRUE)
    if(class(model.inv)=="try-error") {warning("The correlation matrix is singular, an approximation is used")
       ev.mod <- eigen(model)
       ev.mod$values[ev.mod$values < .Machine$double.eps] <- 100 * .Machine$double.eps
       model <- ev.mod$vectors %*% diag(ev.mod$values) %*% t(ev.mod$vectors)
       diag(model)  <- 1
       model.inv <- solve(model)}
    m.inv.r <- model.inv %*% r
    if(is.na(n.obs)) {result$n.obs=NA
    			
    			result$PVAL=NA} else {result$n.obs=n.obs}
    result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n 
    result$criteria <- c("objective"=result$objective,NA,NA)
    if (!is.na(n.obs)) {result$STATISTIC <-  result$objective * ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3)
          if(!is.nan(result$STATISTIC))if (result$STATISTIC <0) {result$STATISTIC <- 0}  
   			if (result$dof > 0) {result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else result$PVAL <- NA}
   	result$Call <- cl
   	
   	#now, find the correlations of the factor scores, even if not estimated, with the factors
   	if(!is.null(phi)) f <- f %*% phi
      w <- try(solve(r,f) ,silent=TRUE)  #these are the factor weights
     if(class(w)=="try-error") {message("In factor.stats, the correlation matrix is singular, an approximation is used")
     ev <- eigen(r)
     ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
       r <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
       diag(r)  <- 1
      w <- solve(r,f)}   #these are the beta weights 
      R2 <- diag(t(w) %*% f)
      r.scores <- cov2cor(t(w) %*% r %*% w)
      result$r.scores <- r.scores 
   	  result$R2 <-R2   #this is the multiple R2 of the scores with the factors
   	  
   	 # result$R2.corrected <- factor.indeterm(r,f)
   	 # result$R2.total <- R2.cor$R2
   	 # result$beta.total <- R2.cor$beta.total
   	  #course coding
   	  keys <- factor2cluster(f) 
   	  covar <- t(keys) %*% r %*% keys 
 
   	  if((nfactors >1) && (dim(covar)[2] >1  )) {
     sd.inv <- diag(1/sqrt(diag(covar)))
     cluster.correl <- sd.inv %*% covar  %*% sd.inv
   	 valid <- t(f) %*% keys %*% sd.inv
   	 result$valid <- diag(valid)
   	 result$score.cor <- cluster.correl} else {sd.inv <- 1/sqrt(covar)
   	                                           if(dim(sd.inv)[1] == 1) sd.inv <- diag(sd.inv)
   	                                           valid <- t(f) %*% keys * sd.inv
   	                                           result$valid <- valid}
   	 result$weights <- w  #the beta weights for factor scores
   	class(result) <- c("psych","stats")
   	return(result)	
   }
   
   
 