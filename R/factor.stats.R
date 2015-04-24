
"factor.stats" <- 
function(r=NULL,f,phi=NULL,n.obs=NA,np.obs=NULL,alpha=.1,fm=NULL) {
   fa.stats(r=r,f=f,phi=phi,n.obs=n.obs,np.obs=np.obs,alpha=alpha,fm=fm)}

"fa.stats" <- 
function(r=NULL,f,phi=NULL,n.obs=NA,np.obs=NULL,alpha=.1,fm=NULL) {
#revised June 21, 2010 to add RMSEA etc. 
#revised August 25, 2011 to add cor.smooth for smoothing
#revised November 10, 2012 to add stats for the minchi option of factoring
#revised February 28, 2014 to emphasize empirical chi 2 and report empirical BIC
#revised March 9, 2015 to report NA if RMSEA values are not in the confidence intervals
cl <- match.call()
conf.level <- alpha 
 if((!is.matrix(f)) && (!is.data.frame(f)))  {#do a number of things that use f as list
   
    if(is.null(r) && (!is.null(f$r)) ) r <- f$r   #we found the correlation while factoring 
    
    #if(is.na(n.obs) && (!is.null(f$np.obs))) {np.obs <- f$np.obs}   
 
 f <- as.matrix(f$loadings)} else {f <- as.matrix(f)}
 
  n <- dim(r)[2]  #number of variables
    if(dim(r)[1] !=n ) {n.obs = dim(r)[1]
                    r <- cor(r,use="pairwise")
                     } 
    if(is.data.frame(r)) r <- as.matrix(r)
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
    if(is.null(np.obs))  {rstar.off <- sum(residual^2)
                          result$ENull <- r2.off * n.obs  #the empirical null model
                          result$chi <- rstar.off * n.obs  #this is the empirical chi square
                          result$rms <- sqrt(rstar.off/(n*(n-1)))  #this is the empirical rmsea
                          result$nh <- n.obs
                             if (result$dof > 0) {result$EPVAL <- pchisq(result$chi, result$dof, lower.tail = FALSE)
                                 result$crms <- sqrt(rstar.off/(2*result$dof) )
                                  result$EBIC <- result$chi - result$dof * log(n.obs) 
                                  result$ESABIC <- result$chi - result$dof * log((n.obs+2)/24) } else {result$EPVAL <- NA
                                 result$crms <- NA
                                  result$EBIC <- NA
                                 result$ESABIC <- NA}
                                             
                          } else {
                           rstar.off <- sum(residual^2 * np.obs)  #weight the residuals by their sample size
                          r2.off <-(r*r * np.obs)   #weight the original by sample size
                          r2.off <- sum(r2.off) -tr(r2.off) 
                          result$chi <- rstar.off  #this is the sample size weighted chi square
                          result$nh <- harmonic.mean(as.vector(np.obs)) #this is the sample weighted cell size
                          result$rms <- sqrt(rstar.off/(result$nh*n*(n-1))) #this is the sample size weighted square root average squared residual
                        if (result$dof > 0) {result$EPVAL <- pchisq(result$chi, result$dof, lower.tail = FALSE)
                                              result$crms <- sqrt(rstar.off/(2*result$nh*result$dof) )
                                              result$EBIC <- result$chi - result$dof * log(result$nh) 
                                              result$ESABIC <- result$chi - result$dof * log((result$nh+2)/24) } else {   #added 2/28/2014 
                                              result$EPVAL <- NA
                                              result$crms <- NA
                                              result$EBIC <- NA
                                              result$ESABIC <- NA
                                              }
                          }

    result$fit <-1-rstar2/r2
    result$fit.off <- 1-rstar.off/r2.off
    result$sd <- sd(as.vector(residual)) #this is the none sample size weighted root mean square residual
    result$factors <- nfactors
    
    result$complexity <- (apply(f,1,function(x) sum(x^2)))^2/apply(f,1,function(x)sum(x^4))
  
    diag(model) <- diag(r)  
    model <- cor.smooth(model)  #this replaces the next few lines with a slightly cleaner approach
    r <- cor.smooth(r)  #this makes sure that the correlation is positive semi-definite
    #although it would seem that the model should always be positive semidefinite so this is probably not necessary
    #cor.smooth approach  added August 25,2011
   #  model.inv <- try(solve(model),silent=TRUE)
   #  if(class(model.inv)=="try-error") {warning("The correlation matrix is singular, an approximation is used")
   #    ev.mod <- eigen(model)
   #   ev.mod$values[ev.mod$values < .Machine$double.eps] <- 100 * .Machine$double.eps
   #   model <- ev.mod$vectors %*% diag(ev.mod$values) %*% t(ev.mod$vectors)
   #    diag(model)  <- 1
   #    #model.inv <- solve(model)
   #    }
   
    m.inv.r <- try(solve(model,r),silent=TRUE) #modified Oct 30, 2009 to perhaps increase precision -- #modified 2015/1/2 to use try
      
     if(class(m.inv.r)=="try-error") {warning("the model inverse times the r matrix is singular, replaced with Identity matrix which means fits are wrong")
            m.inv.r <- diag(1,n,n)}
    if(is.na(n.obs)) {result$n.obs=NA 
    			      result$PVAL=NA} else {result$n.obs=n.obs}
    result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n   #this is what Tucker Lewis call F
    if(is.infinite(result$objective)) {result$objective <- rstar2
                                       message("The determinant of the smoothed correlation was zero.\nThis means the objective function is not defined.\nChi square is based upon observed residuals.")}
    result$criteria <- c("objective"=result$objective,NA,NA)
   
    if (!is.na(n.obs)) {result$STATISTIC <-  chisq <- result$objective * ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3) #from Tucker  and from factanal
    # if (!is.na(n.obs)) {result$STATISTIC <-  chisq <- result$objective * ((n.obs-1)) #from Fox and sem
          if(!is.nan(result$STATISTIC)) if (result$STATISTIC <0) {result$STATISTIC <- 0}  
   			if (result$dof > 0) {result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else {result$PVAL <- NA}
   		}
   	result$Call <- cl
   	
   	#find the Tucker Lewis Index of reliability
   	#Also known as the NNFI which is expressed in terms of Chisq
   	#NNFI <- (chisqNull/dfNull - chisq/df)/(chisqNull/dfNull - 1)
   	#first find the null model 
   	F0 <- sum(diag((r))) - log(det(r)) -n  
   	if(is.infinite(F0))  {F0 <- r2
   	                     message("The determinant of the smoothed correlation was zero.\nThis means the objective function is not defined for the null model either.\nThe Chi square is thus based upon observed correlations.")}
   	Fm <-  result$objective   #objective function of model     
   	Mm <- Fm/( n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2))
   	M0 <- F0* 2 /(n*(n-1))
    nm <- ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3) #
   	result$null.model <- F0
   	result$null.dof <- n * (n-1) /2
   	if (!is.na(n.obs)) {result$null.chisq <-  F0 * ((n.obs-1) -(2 * n + 5)/6 )
                  	result$TLI <- (M0 - Mm)/(M0 - 1/nm)        #NNFI in Fox's sem
                  	if(is.numeric(result$TLI) & !is.nan(result$TLI) & (result$TLI >1)) result$F0 <-1 
    
     #The estimatation of RMSEA and the upper and lower bounds are taken from John Fox's summary.sem with minor modifications
      if(!is.null(result$objective) && (result$dof >0) &&(!is.na(result$objective))) {
      RMSEA <- sqrt(max(result$objective/result$dof - 1/(n.obs-1), 0))        #this is x2/(df*N ) -  1/(N-1)  
   

      
        tail <- conf.level/2    #this had been incorrectly listed as (1-conf.level)/2  which gave extraordinarily narrow confidence boundaries, fixed August 25, 2011
        N <- max <- n.obs
        df <- result$dof
        #chi.sq.statistic <- RMSEA^2 * df * (N - 1) + df
        
        #why isn't this  just chi.sq?
        chi.sq.statistic <- chisq
        
        
        
        
        max <- max(max,chi.sq.statistic) +2* max
        
         #the alternative to this is to use the uniroot technique of Yves Rosseel in  lavaan
         #NOT YET IMPLEMENTED
#        upperlambda <- function(lam)   {tail - pchisq(chi.sq.statistic, df, ncp=lam)^2 }
#       res <- try(uniroot(f = upperlambda,lower=0,upper=max)$root,silent=TRUE) 
#         if(class(res)=="try-error") {message("In factor.stats, I could not find the RMSEA upper bound . Sorry about that")
#                                         res <- NULL
#                                          lam.U <- NA} else {lam.U <- res}
# 			if (is.null(res) || is.na(res$objective) || res$objective < 0){
# 				max <- 0
# 				warning("cannot find upper bound of RMSEA")
# 				break
# 				}				
       
        
        
        while (max > 1){
             res <- try(optimize(function(lam) (tail - pchisq(chi.sq.statistic, df, ncp=lam))^2, interval=c(0, max)),silent=TRUE)
              if(class(res)=="try-error") {message("In factor.stats, I could not find the RMSEA upper bound . Sorry about that")
                                         res <- NULL}
 			if (is.null(res) || is.na(res$objective) || res$objective < 0){
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
       if(is.na(max)) max <- N
         while (max > 1){#        this just iterates in to get a value
             res <- try(optimize(function(lam) (1 - tail - pchisq(chi.sq.statistic, df, ncp=lam))^2, interval=c(0, max)),silent=TRUE)
              if(class(res)=="try-error") {message("In factor.stats, I could not find the RMSEA lower bound. Sorry about that")
                                         res <- NULL}
             if (is.null(res)) {break}
             if (sqrt(res$objective) < tail/100) break
             max <- max/2
 			if (is.na(res$objective) || res$objective < 0){
 				max <- 0
 				warning("cannot find lower bound of RMSEA")
 				break
 				}				
             }


        lam.L <- if (max <= 1) NA else res$minimum  #lam is the ncp
       #this RMSEA calculation is probably not right because it will sometimes (but rarely) give cis that don't include the estimate   
        RMSEA.U <- sqrt(lam.U/((N)*df) )   #lavaan uses sqrt(lam.U/((N)*df) )  sem uses sqrt(lam.U/((N-1)*df) )
        RMSEA.L <- min(sqrt(lam.L/((N)*df) ),RMSEA)
        if(!is.na(RMSEA.U) && RMSEA.U < RMSEA) RMSEA.U <- NA 
          if(!is.na(RMSEA.L) && RMSEA.L  > RMSEA) RMSEA.L  <- NA 
       result$RMSEA <- c(RMSEA, RMSEA.L, RMSEA.U, conf.level)
       names(result$RMSEA) <- c("RMSEA","lower","upper","confidence")
       result$BIC <- chisq - df * log(N) 
       result$SABIC <- chisq - df * log((N+2)/24)  #added 1/27/2014
        }
      }  
  
        
   	
   	#now, find the correlations of the factor scores, even if not estimated, with the factors
   	if(!is.null(phi)) f <- f %*% phi   #convert the pattern to structure coefficients
   	 r <- cor.smooth(r)
      w <- try(solve(r,f) ,silent=TRUE)  #these are the regression factor weights
     if(class(w)=="try-error") {message("In factor.stats, the correlation matrix is singular, an approximation is used")
     ev <- eigen(r)
     if(is.complex(ev$values)) {warning("complex eigen values detected by factor stats, results are suspect")
                
                 } else { 
     ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
       r <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
       diag(r)  <- 1
     w <- try(solve(r,f) ,silent=TRUE)  #these are the factor weights
     if(class(w)=="try-error") {warning("In factor.stats, the correlation matrix is singular, and we could not calculate the beta weights for factor score estimates")
     w <- diag(1,dim(r)[1])
     }   #these are the beta weights 
    }}
      R2 <- diag(t(w) %*% f)
      if(is.null(fm)) {
     if(prod(R2) < 0 ) {message("In factor.stats: The factor scoring weights matrix is probably singular -- Factor score estimate results are likely incorrect.\n Try a different factor extraction method\n")
                      R2[abs(R2) > 1] <- NA
                      R2[R2 <= 0] <- NA
                     }
     if ((max(R2,na.rm=TRUE) > (1 + .Machine$double.eps)) ) {message("The estimated weights for the factor scores are probably incorrect.  Try a different factor extraction method.")}
     }
      r.scores <- cov2cor(t(w) %*% r %*% w) 
      result$r.scores <- r.scores 
   	  result$R2 <- R2   #this is the multiple R2 of the scores with the factors
   	  
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
   	                                           valid <- try(t(f) %*% keys * sd.inv)
   	                                           result$valid <- valid}
   	 result$weights <- w  #the beta weights for factor scores
   	class(result) <- c("psych","stats")
   	return(result)	
   }
   
   
 