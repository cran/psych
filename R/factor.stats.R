"factor.stats" <- 
function(r,f,phi=NULL,n.obs=NA,alpha=.1) {
#revised June 21, 2010 to add RMSEA etc. 

cl <- match.call()
conf.level <- alpha 
 n <- dim(r)[2]  #number of variables
 if(dim(r)[1] !=n ) {n.obs = dim(r)[1]
                    r <- cor(r,use="pairwise")
                     }
 if(is.data.frame(r)) r <- as.matrix(r)
 if((!is.matrix(f)) && (!is.data.frame(f)))  {f <- as.matrix(f$loadings)} else {f <- as.matrix(f)}
 nfactors <- dim(f)[2]  # number of factors
 if(is.null(phi)) {model <- f %*%  t(f)} else {model <- f %*% phi %*% t(f)}

    residual<- r - model
  
    r2 <- sum(r*r)
    rstar2 <- sum(residual*residual)
    result <- list(residual = residual)
 
    result$dof <- dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
    #r2.off <- r
    #diag(r2.off) <- 0
   # r2.off <- sum(r2.off^2)
   r2.off <- r2 - tr(r)
    diag(residual) <- 0
    rstar.off <- sum(residual^2)
    result$fit <-1-rstar2/r2
    result$fit.off <- 1-rstar.off/r2.off
    result$sd <- sd(as.vector(residual))
    if(dof > 0) result$crms <- sqrt(rstar.off/(2*dof))  #this is the empirical rmsea
    if(dof > 0) result$rms <- sqrt(rstar.off/(2*n*(n-1)))  #this is the empirical rmsea
    
    
    
    result$factors <- nfactors
  
    diag(model) <- diag(r)   
     model.inv <- try(solve(model),silent=TRUE)
    if(class(model.inv)=="try-error") {warning("The correlation matrix is singular, an approximation is used")
       ev.mod <- eigen(model)
       ev.mod$values[ev.mod$values < .Machine$double.eps] <- 100 * .Machine$double.eps
       model <- ev.mod$vectors %*% diag(ev.mod$values) %*% t(ev.mod$vectors)
       diag(model)  <- 1
       #model.inv <- solve(model)
       }
   
     m.inv.r <- solve(model,r)  #modified Oct 30, 2009 to perhaps increase precision
    if(is.na(n.obs)) {result$n.obs=NA 
    			      result$PVAL=NA} else {result$n.obs=n.obs}
    result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n   #this is what Tucker Lewis call F
    result$criteria <- c("objective"=result$objective,NA,NA)
   
    if (!is.na(n.obs)) {result$STATISTIC <-  chisq <- result$objective * ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3) #from Tucker  and from factanal
    # if (!is.na(n.obs)) {result$STATISTIC <-  chisq <- result$objective * ((n.obs-1)) #from Fox and sem
          if(!is.nan(result$STATISTIC)) if (result$STATISTIC <0) {result$STATISTIC <- 0}  
   			if (result$dof > 0) {result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else result$PVAL <- NA}
   	result$Call <- cl
   	
   	#find the Tucker Lewis Index of reliability
   	#Also known as the NNFI which is expressed in terms of Chisq
   	#NNFI <- (chisqNull/dfNull - chisq/df)/(chisqNull/dfNull - 1)
   	#first find the null model 
   	F0 <- sum(diag((r))) - log(det(r)) -n     
   	Fm <-  result$objective   #objective function of model     
   	Mm <- Fm/( n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2))
   	M0 <- F0* 2 /(n*(n-1))
    nm <- ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3) #
   	result$null.model <- F0
   	result$null.dof <- n * (n-1) /2
   	if (!is.na(n.obs)) {result$null.chisq <-  F0 * ((n.obs-1) -(2 * n + 5)/6 )
                  	result$TLI <- (M0 - Mm)/(M0 - 1/nm)        #NNFI in Fox's sem
                  	if(is.numeric(result$TLI) & !is.nan(result$TLI) & (result$TLI >1)) result$F0 <-1 
     
     #The estimatation of RMSEA and the upper and lower bounds is taken from John Fox's sem with minor modifications
      if(!is.null(result$objective) && (result$dof >0)) {
      RMSEA <- sqrt(max(result$objective/result$dof - 1/(n.obs-1), 0))        #this is x2/(df * (N-1)  - 1/(N-1)
   

      
        tail <- (1 - conf.level)/2 
        N <- max <- n.obs
        df <- result$dof
        chi.sq.statistic <- RMSEA^2 * df * (N - 1) + df
        max <- max(max,chi.sq.statistic)
        while (max > 1){
            res <- optimize(function(lam) (tail - pchisq(chi.sq.statistic, df, ncp=lam))^2, interval=c(0, max))
			if (is.na(res$objective) || res$objective < 0){
				max <- 0
				warning("cannot find upper bound of RMSEA")
				break
				}				
            if (sqrt(res$objective) < tail/100) break
            max <- max/2
            }
        lam.U <- if (max <= 1) NA else res$minimum
      # max <- max(max,lam.U)
      max <- lam.U
        while (max > 1){
            res <- optimize(function(lam) (1 - tail - pchisq(chi.sq.statistic, df, ncp=lam))^2, interval=c(0, max))
            if (sqrt(res$objective) < tail/100) break
            max <- max/2
			if (is.na(res$objective) || res$objective < 0){
				max <- 0
				warning("cannot find lower bound of RMSEA")
				break
				}				
            }
        lam.L <- if (max <= 1) NA else res$minimum  #lam is the ncp
       #this RMSEA calculation is not right   
        RMSEA.U <- sqrt(lam.U/((N-1)*df) )
        RMSEA.L <- min(sqrt(lam.L/((N-1)*df) ),RMSEA)
       result$RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
       names(result$RMSEA) <- c("RMSEA","lower","upper","confidence")
        result$BIC <- chisq - df * log(N) 
        }
      }  
       
        
   	
   	#now, find the correlations of the factor scores, even if not estimated, with the factors
   	if(!is.null(phi)) f <- f %*% phi   #convert the pattern to structure coefficients
      w <- try(solve(r,f) ,silent=TRUE)  #these are the factor weights
     if(class(w)=="try-error") {message("In factor.stats, the correlation matrix is singular, an approximation is used")
     ev <- eigen(r)
     ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
       r <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
       diag(r)  <- 1
     w <- try(solve(r,f) ,silent=TRUE)  #these are the factor weights
     if(class(w)=="try-error") {warning("In factor.stats, the correlation matrix is singular, and we could not calculate the beta weights for factor score estimates")
     w <- diag(1,dim(r)[1])
     }   #these are the beta weights 
    }
      R2 <- diag(t(w) %*% f)
     if(prod(R2) <0 ) {message("The matrix is probably singular -- Factor score estimate results are likely incorrect")
                      R2[abs(R2) > 1] <- NA
                      #added to 
                     }
     #if ((max(R2) > (1 + .Machine$double.eps)) ) {message("The estimated weights for the factor scores are probably incorrect.  Try a different factor extraction method.")}
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
   
   
 