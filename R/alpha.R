"alpha" <- 
    function(x,keys=NULL,cumulative=FALSE,title=NULL,max=10,na.rm=TRUE,check.keys=FALSE,n.iter=1,delete=TRUE,use="pairwise",warnings=TRUE,n.obs=NULL) {  #find coefficient alpha given a data frame or a matrix
    
    alpha.1 <- function(C,R) {
    n <- dim(C)[2]
    alpha.raw <- (1- tr(C)/sum(C))*(n/(n-1))
    sumR <- sum(R)
    alpha.std <-  (1- n/sumR)*(n/(n-1))
    smc.R <- smc(R)
    G6 <- (1- (n-sum(smc.R))/sumR)
    av.r <- (sumR-n)/(n*(n-1))
    mod1 <- matrix(av.r,n,n)
    Res1 <- R - mod1
    GF1 =  1- sum(Res1^2)/sum(R^2)
    Rd <- R - diag(R)
    diag(Res1) <- 0
    GF1.off <- 1 - sum(Res1^2)/sum(Rd^2)  
    sn <- n*av.r/(1-av.r)
   # Q = (2 * n^2/((n-1)^2*(sum(C)^3))) * (sum(C) * (tr(C^2) + (tr(C))^2) - 2*(tr(C) * sum(C^2))) #corrected 1/15/16 
    Q = (2 * n^2/((n - 1)^2 * (sum(C)^3))) * (sum(C) * (tr(C%*%C) +  (tr(C))^2) - 2 * (tr(C) * sum(C%*%C)))   #correction from Tamaki Hattori
    result <- list(raw=alpha.raw,std=alpha.std,G6=G6,av.r=av.r,sn=sn,Q=Q,GF1,GF1.off)
    return(result)
    }
    
    #begin main function
    cl <- match.call()
    if(!is.matrix(x) && !is.data.frame(x)) stop('Data must either be a data frame or a matrix')
    nvar <- dim(x)[2]
    nsub <- dim(x)[1]
    scores <- NULL
    response.freq <- NULL
   
    if (nsub !=nvar)  { #find the correlations if we are given  raw data
       item.var <- apply(x,2,sd,na.rm=na.rm)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0) && delete) {
            for (baddy in 1:length(bad)) {warning( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
            nvar <- nvar - length(bad)
             }
         response.freq <- response.frequencies(x,max=max)
         C <- cov(x,use=use)} else {C <- x}
        
        if(is.null(colnames(x)))  colnames(x) <- paste0("V",1:nvar)
       
         #if(check.keys && is.null(keys)) {
            p1 <- principal(x,scores=FALSE)
               if(any(p1$loadings < 0)) {if (check.keys) {if(warnings) warning("Some items were negatively correlated with total scale and were automatically reversed.\n This is indicated by a negative sign for the variable name.") 
                    keys <- 1- 2* (p1$loadings < 0 ) } else {
                       if(is.null(keys) && warnings ) {warning("Some items were negatively correlated with the total scale and probably \nshould be reversed.  \nTo do this, run the function again with the 'check.keys=TRUE' option")
                       if(warnings) cat("Some items (",rownames(p1$loadings)[(p1$loadings < 0)],") were negatively correlated with the total scale and \nprobably should be reversed.  \nTo do this, run the function again with the 'check.keys=TRUE' option")
                        }} 
            }  #keys is now a vector of 1s and -1s
           # }
            
         if(is.null(keys)) {keys <- rep(1,nvar)} else {  
         			keys<- as.vector(keys) 
         			if(length(keys) < nvar) {temp <- keys  #this is the option of keying just the reversals
         			                         keys <- rep(1,nvar)
         			                         names(keys) <- colnames(x)
         			                         keys[temp] <- -1}
         			                   } 
        			 key.d <- diag(keys)
                     C <- key.d %*% C %*% key.d
                     signkey <- strtrim(keys,1)
            		 signkey[signkey=="1"] <- ""
                     colnames(x) <- paste(colnames(x),signkey,sep="")
                     
      if (nsub !=nvar)  {   #raw data      
         	if(any(keys < 0 )) { 
         	     min.item <- min(x,na.rm=na.rm)
                 max.item <- max(x,na.rm=na.rm)
         	   adjust <- max.item + min.item
         	   flip.these <- which(keys < 0 )
         	   x[,flip.these]  <- adjust - x[,flip.these] 
         	    }
         	

        if(cumulative) {total <- rowSums(x,na.rm=na.rm) } else {total <- rowMeans(x,na.rm=na.rm)}
                mean.t <- mean(total,na.rm=na.rm)
                sdev <- sd(total,na.rm=na.rm) 
                raw.r <- cor(total,x,use=use)

                       
        t.valid <- colSums(!is.na(x))} else {   #we are working with a correlation matrix
                 total <- NULL
                 totals <- TRUE
                 }
          
         R <- cov2cor(C)
         drop.item <- vector("list",nvar)
         alpha.total <- alpha.1(C,R)
         if(nvar > 2) {
         for (i in 1:nvar) {
         drop.item[[i]] <- alpha.1(C[-i,-i,drop=FALSE],R[-i,-i,drop=FALSE])
                            } 
         } else {drop.item[[1]] <- drop.item[[2]] <- c(rep(R[1,2],2),smc(R)[1],R[1,2],NA,NA)}
        by.item <- data.frame(matrix(unlist(drop.item),ncol=8,byrow=TRUE))
        
                  #allows us to specify the number of subjects for correlation matrices
        if(max(nsub,n.obs) > nvar) {by.item[6] <- sqrt(by.item[6]/(max(nsub,n.obs)) )
         by.item <- by.item[-c(7:8)]
         colnames(by.item) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","alpha se") } else {
         
             by.item <- by.item[-c(6:8)]
             colnames(by.item) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N") }
        rownames(by.item) <- colnames(x)
        
        Vt <- sum(R)
        item.r <- colSums(R)/sqrt(Vt)  #this is standardized r
       
     #correct for item overlap by using  smc 
        RC <-R
        diag(RC) <-smc(R)
        Vtc <- sum(RC)
        item.rc <-colSums(RC)/sqrt(Vtc)
     #yet one more way to correct is to correlate item with rest of scale
      if(nvar > 1) {
      r.drop <- rep(0,nvar)
        for (i in 1:nvar) { v.drop <- sum(C[-i,-i,drop=FALSE])
          c.drop <- sum(C[,i]) - C[i,i]
          r.drop[i] <- c.drop/sqrt(C[i,i]*v.drop)
              }
      }
     
     #  
        item.means <- colMeans(x, na.rm=na.rm )
        item.sd <-  apply(x,2,sd,na.rm=na.rm)
        if(nsub > nvar) {
         Unidim <- alpha.total[7]
        	 Fit.off <- alpha.total[8]
            ase = sqrt(alpha.total$Q/nsub)
        	alpha.total <- data.frame(alpha.total[1:5],ase=ase,mean=mean.t,sd=sdev)
           colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","ase","mean","sd")
        	
        	  
        	 alpha.total <- data.frame(alpha.total[1:5],ase=ase,mean=mean.t,sd=sdev)
        	colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","ase","mean","sd")
        	rownames(alpha.total) <- ""
        	stats <- data.frame(n=t.valid,raw.r=t(raw.r),std.r =item.r,r.cor = item.rc,r.drop = r.drop,mean=item.means,sd=item.sd)
        	} else {
        	if(is.null(n.obs)) {
        	       Unidim <- alpha.total[7]
        	       Fit.off <- alpha.total[8]
        	      alpha.total <- data.frame(alpha.total[1:5])  #fixed 27/7/14 
        	        colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)" ,"average_r","S/N") } else {
        	        Unidim <- alpha.total[7]
        	 		Fit.off <- alpha.total[8] 
        	        alpha.total <- data.frame(alpha.total[1:5],ase=sqrt(alpha.total$Q/n.obs))
        	        colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)" ,"average_r","S/N","ase")}
        	        rownames(alpha.total) <- "" 
        	        

                  stats <- data.frame(r =item.r,r.cor = item.rc,r.drop = r.drop) #added r.drop 10/12/13
        	}
       	rownames(stats) <- colnames(x)
       	
       	#added measures of unidimensionality  Feb 24, 2016
        #found in alpha.1
        # the basic idea is how big are the residuals given the model
        #we can compare them to the total R or the off diagonal R.
      
       	#end of unidimensionality statistics
       	if(n.iter > 1) {#do a bootstrap confidence interval for alpha
   #    	 if(!require(parallel)) {message("The parallel package needs to be installed to run mclapply")}
        if(nsub == nvar) {message("bootstrapped confidence intervals require raw data") 
                          boot <- NULL
                          boot.ci <- NULL } else {
       	 boot <- vector("list",n.iter)
       	 boot <- mclapply(1:n.iter,function(XX) {
       	 xi <- x[sample.int(nsub,replace=TRUE),]
       	   C <- cov(xi,use="pairwise")
       	           if(!is.null(keys)) {key.d <- diag(keys)
                                      xi <- key.d %*% C %*% key.d}
                
                                      
         R <- cov2cor(C)
       	 alpha.1(C,R)
       	 })  #end of mclapply 
       
       	  boot <- matrix(unlist(boot),ncol=8,byrow=TRUE)
       	  colnames(boot) <- c("raw_alpha","std.alpha","G6(smc)","average_r","s/n","ase","Unidim","Goodfit")
       	  boot.ci <- quantile(boot[,1],c(.025,.5,.975))
       	}} else {boot=NULL
       	         boot.ci <- NULL}
       	names(Unidim) <- "Unidim"
       	names(Fit.off) <- "Fit.off" 
        result <- list(total=alpha.total,alpha.drop=by.item,item.stats=stats,response.freq=response.freq,keys=keys,scores = total,nvar=nvar,boot.ci=boot.ci,boot=boot,Unidim=Unidim,Fit=Fit.off,call=cl,title=title)
        class(result) <- c("psych","alpha")
        return(result) 
    }
  #modified Sept 8, 2010 to add r.drop feature  
  #modified October 12, 2011 to add apply to the sd function
  #modified November 2, 2010 to use sd instead of SD
  #January 30, 2011  - added the max category parameter (max)
  #June 20, 2011 -- revised to add the check.keys option as suggested by Jeremy Miles
  #Oct 3, 2013 check for variables with no variance and drop them with a warning
  #November 22, 2013  Added the standard error as suggested by 
  #modified December 6, 2013 to add empirical confidence estimates
  #modified January 9, 2014 to add multicore capabilities to the bootstrap 
  #corrected December 18 to allow reverse keying for correlation matrices as well as raw data
  #modified 1/16/14 to add S/N to summary stats
  #added item.c  (raw correlation) 1/10/15
  #corrected 1/16/16 corrected the formula for Q following a suggestion by Tamaki Hattori
  #added the n.obs option to allow us to find standard errors even from correlation matrices