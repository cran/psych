"alpha" <- 
    function(x,keys=NULL,cumulative=FALSE,title=NULL,max=10,na.rm=TRUE,check.keys=FALSE,n.iter=1,delete=TRUE,use="pairwise",warnings=TRUE,n.obs=NULL,impute=NULL) {  #find coefficient alpha given a data frame or a matrix
    
 alpha.1 <- function(C,R) {
    n <- dim(C)[2]
    alpha.raw <- (1- tr(C)/sum(C))*(n/(n-1))
    sumR <- sum(R)
    alpha.std <-  (1- n/sumR)*(n/(n-1))  
    smc.R <- smc(R)
    G6 <- (1- (n-sum(smc.R))/sumR)
    av.r <- (sumR-n)/(n*(n-1))
    R.adj <- R[lower.tri(R)]
   # diag(R.adj) <- NA
   # var.r  <- var(as.vector(R.adj),na.rm=TRUE)
   var.r <- var(R.adj,na.rm=TRUE)  #
   med.r <- median(R.adj,na.rm=TRUE)  #added 4/22/18
    mod1 <- matrix(av.r,n,n)
    Res1 <- R - mod1
    GF1 =  1- sum(Res1^2)/sum(R^2)
    Rd <- R - diag(R)
    diag(Res1) <- 0
    GF1.off <- 1 - sum(Res1^2)/sum(Rd^2)  
    sn <- n*av.r/(1-av.r)
   # Q = (2 * n^2/((n-1)^2*(sum(C)^3))) * (sum(C) * (tr(C^2) + (tr(C))^2) - 2*(tr(C) * sum(C^2))) #corrected 1/15/16 
    Q = (2 * n^2/((n - 1)^2 * (sum(C)^3))) * (sum(C) * (tr(C%*%C) +  (tr(C))^2) - 2 * (tr(C) * sum(C%*%C)))   #correction from Tamaki Hattori
    result <- list(raw=alpha.raw,std=alpha.std,G6=G6,av.r=av.r,sn=sn,Q=Q,GF1,GF1.off,var.r = var.r,med.r=med.r)
    return(result)
    }
    
#begin main function
    cl <- match.call()
    if(!is.matrix(x) && !is.data.frame(x)) stop('Data must either be a data frame or a matrix')
    if(!inherits(x[1], "data.frame")) x <- fix.dplyr(x)    #to get around a problem created by dplyr
    if(!is.null(keys)){# 4 cases  1 it is a list, 2 is a vector of character, 3 it is keys matrix  4 it is a list of items to  reverse
    if( is.list(keys)) { select <- sub("-","",unlist(keys))   #added 9/26/16 to speed up scoring one scale from many
      	x <- x[,select] 
      	keys <- make.keys(x,keys)} else {if(!is.numeric(keys)){
      	temp <- rep(1,ncol(x))
      	temp[(colnames(x) %in% keys)] <- -1  #note that if the keys vect has negatives, they will not be reversed, but rather the positive ones will, net result is the same, but item correlations are backward
      	keys <- temp
      }}
      }
    if(missing(impute) | (!is.null(impute) && (impute =="none"))) impute <- NULL  
    nvar <- dim(x)[2]
    nsub <- dim(x)[1]
    scores <- NULL
    response.freq <- NULL
    raw <- FALSE
    if (!isCovariance(x))  { #find the correlations if we are given  raw data
      raw <- TRUE
       if(!is.null(impute) ) {
       		if(impute=="median") {item.impute <- apply(x,2,median,na.rm=na.rm)} else {item.impute <- apply(x,2,mean,na.rm=na.rm) } #column values
 
            find.na <- which(is.na(x),arr.ind=TRUE)
            if(NROW(find.na) >0) x[find.na] <- item.impute[find.na[,2]]   #added 08/28/20 to fix a problem where we putting in imputed values for all the variables
            }
       item.var <- apply(x,2,sd,na.rm=na.rm)
       bad <- which((item.var <= 0)|is.na(item.var))
       n.bad <-length(bad)
       if((n.bad > 0) && delete) {
            for (baddy in 1:n.bad) {warning( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted but still is counted in the score")}
            first.bad <- x[1,bad]
            x <- x[,-bad] 
            nvar <- nvar - n.bad
             }
         response.freq <- response.frequencies(x,max=max)
         C <- cov(x,use=use)} else {C <- x}  #to find average R, we need to find R, not C  cov2cor(C) != R  !
        
        if(is.null(colnames(x)))  colnames(x) <- paste0("V",1:nvar)
         #flip items if needed and wanted
         #if(check.keys && is.null(keys)) {
            p1 <- principal(x,scores=FALSE)

               if(any(p1$loadings < 0)) {if (check.keys) {if(warnings) warning("Some items were negatively correlated with the first principal component and were automatically reversed.\n This is indicated by a negative sign for the variable name.") 
                    keys <- 1- 2* (p1$loadings < 0)
                      } else {
                       if(is.null(keys) && warnings ) {warning("Some items were negatively correlated with the first principal component and probably \nshould be reversed.  \nTo do this, run the function again with the 'check.keys=TRUE' option")
                       if(warnings) cat("Some items (",rownames(p1$loadings)[(p1$loadings < 0)],") were negatively correlated with the first principal component and \nprobably should be reversed.  \nTo do this, run the function again with the 'check.keys=TRUE' option")
                        keys <- rep(1,nvar)
                         } 
                   
                    }
            }  #keys is now a vector of 1s and -1s
           #names(keys) <- colnames(x)
           # }
            
         if(is.null(keys)) {keys <- rep(1,nvar)
              
                 names(keys) <- colnames(x)} else {  
         			keys<- as.vector(keys) 
         	        names(keys) <- colnames(x)
         			if(length(keys) < nvar) {temp <- keys  #this is the option of keying just the reversals
         			                         keys <- rep(1,nvar)
         			                         names(keys) <- colnames(x)
         			                         keys[temp] <- -1
         			                      }
         			                         
         			                   } 
         			           
         			           
        			 key.d <- diag(keys)
                     C <- key.d %*% C %*% key.d
                     signkey <- strtrim(keys,1)
            		 signkey[signkey=="1"] <- ""
                     colnames(x) <- paste(colnames(x),signkey,sep="")
                     
      if (raw)  {   #raw data      
         	if(any(keys < 0 )) { 
         	     min.item <- min(x,na.rm=na.rm)
                 max.item <- max(x,na.rm=na.rm)
         	   adjust <- max.item + min.item
         	   flip.these <- which(keys < 0 )
         	   x[,flip.these]  <- adjust - x[,flip.these] 
         	  #note that we do not reverse key bad items 
         	    }
         	
        if(n.bad > 0 ){sum.bad <- sum(first.bad, na.rm=TRUE) 
                      if(cumulative) {total <- rowSums(x,na.rm=na.rm) + sum.bad  } else {total <- (nvar * rowMeans(x,na.rm=na.rm) + sum.bad) /(nvar+n.bad)}} else {
        if(cumulative) {total <- rowSums(x,na.rm=na.rm) } else {total <- rowMeans(x,na.rm=na.rm)}
        }
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
         } else {drop.item[[1]]<- c(C[1,2]/C[2,2],R[1,2],smc(R)[1],R[1,2],R[1,2]/(1-R[1,2]),NA,NA,NA,0,R[1,2])
                 drop.item[[2]] <- c(C[1,2]/C[1,1],R[1,2],smc(R)[1],R[1,2],R[1,2]/(1-R[1,2]),NA,NA,NA,0,R[1,2]) }  #added the extra 2 NA June 18, 2017  -- and further fixed May 7, 2020
       
        by.item <- data.frame(matrix(unlist(drop.item),ncol=10,byrow=TRUE)) 
        
                  #allows us to specify the number of subjects for correlation matrices
        if(max(nsub,n.obs) > nvar) {by.item[6] <- sqrt(by.item[6]/(max(nsub,n.obs)) )
         by.item <- by.item[-c(7:8)]
         colnames(by.item) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","alpha se","var.r","med.r") } else {
         
             by.item <- by.item[-c(6:8)]
             colnames(by.item) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","var.r","med.r") }
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
        if(raw) {
         Unidim <- alpha.total[7]
         var.r <- alpha.total[[9]]
        	 Fit.off <- alpha.total[8]
            ase = sqrt(alpha.total$Q/nsub)
        	#alpha.total <- data.frame(alpha.total[1:5],ase=ase,mean=mean.t,sd=sdev)
          # colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","ase","mean","sd")
        	
        	  
        	 alpha.total <- data.frame(alpha.total[1:5],ase=ase,mean=mean.t,sd=sdev,med.r =alpha.total[10])
        	colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)","average_r","S/N","ase","mean","sd","median_r")
        	rownames(alpha.total) <- ""
        	stats <- data.frame(n=t.valid,raw.r=t(raw.r),std.r =item.r,r.cor = item.rc,r.drop = r.drop,mean=item.means,sd=item.sd)
        	} else {
        	if(is.null(n.obs)) {
        	       Unidim <- alpha.total[7]
        	       Fit.off <- alpha.total[8]
        	       var.r <- alpha.total[9] 
        	       med.r <- alpha.total[10]
        	      alpha.total <- data.frame(alpha.total[c(1:5,10)])  #fixed 27/7/14 
        	        colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)" ,"average_r","S/N","median_r") } else {
        	        Unidim <- alpha.total[7]
        	 		Fit.off <- alpha.total[8] 
        	 		        	       var.r <- alpha.total[9] 
        	        alpha.total <- data.frame(alpha.total[1:5],ase=sqrt(alpha.total$Q/n.obs),alpha.total[10])
        	        colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)" ,"average_r","S/N","ase","median_r")}
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
        if(!raw) {message("bootstrapped confidence intervals require raw data") 
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
       
       	  boot <- matrix(unlist(boot),ncol=10,byrow=TRUE)
       	  colnames(boot) <- c("raw_alpha","std.alpha","G6(smc)","average_r","s/n","ase","Unidim","Goodfit","var.r","median.r")
       	  boot.ci <- quantile(boot[,1],c(.025,.5,.975))
       	}} else {boot=NULL
       	         boot.ci <- NULL}
       	names(Unidim) <- "Unidim"
       	names(Fit.off) <- "Fit.off" 
       	feldt <- alpha.ci(alpha.total[1],nsub,nvar)
       	keys <- keys2list(as.matrix(keys))
        result <- list(total=alpha.total,alpha.drop=by.item,item.stats=stats,response.freq=response.freq,keys=keys,scores = total,nvar=nvar,boot.ci=boot.ci,boot=boot,feldt=feldt,Unidim=Unidim,var.r=var.r,Fit=Fit.off,call=cl,title=title)
        class(result) <- c("psych","alpha")
        return(result) 
    }
  #modified Sept 8, 2010 to add r.drop feature  
  #modified October 12, 2011 to add apply to the sd function
  #modified November 2, 2010 to use sd instead of SD
  #January 30, 2011  - added the max category parameter (max)
  #June 20, 2011 -- revised to add the check.keys option as suggested by Jeremy Miles
  #Oct 3, 2013 check for variables with no variance and drop them with a warning
  #May 27, 2022 convert the keys output to a named list so that we can use if in other functions
  #November 22, 2013  Added the standard error as suggested by 
  #modified December 6, 2013 to add empirical confidence estimates
  #modified January 9, 2014 to add multicore capabilities to the bootstrap 
  #corrected December 18 to allow reverse keying for correlation matrices as well as raw data
  #modified 1/16/14 to add S/N to summary stats
  #added item.c  (raw correlation) 1/10/15
  #corrected 1/16/16 corrected the formula for Q following a suggestion by Tamaki Hattori
  #added the n.obs option to allow us to find standard errors even from correlation matrices
  #7/4/23 changed the reference to total score to be first principal component
  
#a kludge to get around a problem introduced by dplyr which changes the class structure of data frames.
#created in response to a problem raised by Adam Liter   (February, 2017) 
#probably not necessary anymore if we use inherits(x,"data.frame")
"fix.dplyr" <- function (object) {
   if(is.data.frame(object)) {   
   cn <- class(object)
df <- which(cn=="data.frame")
cn.not <- cn[-df]
cn <- c("data.frame",cn.not)
class(object) <- cn
      } 
invisible(object) 
}

#apply the Duhacheck and Iacobucci estimates
#compare with Feldt's estimate
"alpha.ci" <- function(alpha,n.obs,n.var=NULL,p.val=.05,digits=2) {
# Q = (2 * n^2/((n - 1)^2 * (sum(C)^3))) * (sum(C) * (tr(C%*%C) +  (tr(C))^2) - 2 * (tr(C) * sum(C%*%C)))   #correction from Tamaki Hattori
 CI.high <- 1- (1-alpha)* qf(p.val/2,df1=n.obs-1,df2=(n.obs-1)*(n.var-1))
 CI.low <- 1-  (1-alpha)* qf(1-p.val/2,df1=n.obs-1,df2=(n.obs-1)*(n.var-1))  #Inf changed inf to nvar -1  3/3/22
#CI.high <- 1- (1-alpha)* qf(p.val/2,df1=n.obs-1,df2=Inf)
#CI.low <- 1-  (1-alpha)* qf(1-p.val/2,df1=n.obs-1,df2=Inf)
if(!is.null(n.var)) {r.bar <- alpha/(n.var - alpha*(n.var-1)) } else {r.bar=NA}
result <- list(lower.ci =CI.low,alpha=alpha,upper.ci=CI.high,r.bar=r.bar)
#print(result,digits=digits)
class(result) <- c("psych","alpha.ci")
invisible(result)
}

#convert alpha the the average r
alpha2r <- function(alpha,n.var) {
   r.bar <- alpha/(n.var - alpha*(n.var-1))
  r.bar}