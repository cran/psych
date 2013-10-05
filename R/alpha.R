"alpha" <- 
    function(x,keys=NULL,cumulative=FALSE,title=NULL,max=10,na.rm=TRUE,check.keys=TRUE,delete=TRUE) {  #find coefficient alpha given a data frame or a matrix
    
    alpha.1 <- function(C,R) {
    n <- dim(C)[2]
    alpha.raw <- (1- tr(C)/sum(C))*(n/(n-1))
    alpha.std <-  (1- n/sum(R))*(n/(n-1))
    smc.R <- smc(R)
    G6 <- (1- (n-sum(smc.R))/sum(R))
    av.r <- (sum(R)-n)/(n*(n-1))
    result <- list(raw=alpha.raw,std=alpha.std,G6,av.r)
    return(result)
    }
    
    #begin main function
    cl <- match.call()
    if(!is.matrix(x) && !is.data.frame(x)) stop('Data must either be a data frame or a matrix')
    nvar <- dim(x)[2]
    nsub <- dim(x)[1]
    scores <- NULL
    response.freq <- NULL
   
   
    if (nsub !=nvar)  {
       item.var <- apply(x,2,sd,na.rm=na.rm)
       bad <- which((item.var==0)|is.na(item.var))
       if((length(bad) > 0) && delete) {
            for (baddy in 1:length(bad)) {warning( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
            nvar <- nvar - length(bad)
             }
         response.freq <- response.frequencies(x,max=max)
         C <- cov(x,use="pairwise") 
         
         if(check.keys && is.null(keys)) {
            p1 <- principal(x)
            if(any(p1$loadings < 0)) warning("Some items were negatively correlated with total scale and were automatically reversed.")
            keys <- 1- 2* (p1$loadings < 0 ) 
            }  #keys is now a vector of 1s and -1s
            
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
            min.item <- min(x,na.rm=na.rm)
         	max.item <- max(x,na.rm=na.rm)
         	if(any(keys <0 )) { 
         	   adjust <- max.item + min.item
         	   flip.these <- which(keys < 0 )
         	   x[,flip.these]  <- adjust - x[,flip.these] 
         	    }
         	

        if(cumulative) {total <- rowSums(x,na.rm=na.rm) } else {total <- rowMeans(x,na.rm=na.rm)}
                   mean.t <- mean(total,na.rm=na.rm)
                  sdev <- sd(total,na.rm=na.rm) 
                       
        t.valid <- colSums(!is.na(x))} else {   #we are working with a correlation matrix
                 total <- NULL
                 totals <- TRUE
                 C <- as.matrix(x) 
                  if(!is.null(keys)) {key.d <- diag(keys)
                                      C <- key.d %*% C %*% key.d}}
         R <- cov2cor(C)
         drop.item <- list()
         alpha.total <- alpha.1(C,R)
         if(nvar > 2) {
         for (i in 1:nvar) {
         drop.item[[i]] <- alpha.1(C[-i,-i,drop=FALSE],R[-i,-i,drop=FALSE])
                            } 
         } else {drop.item[[1]] <- drop.item[[2]] <- c(rep(R[1,2],2),smc(R)[1],R[1,2])}
        by.item <- data.frame(matrix(unlist(drop.item),ncol=4,byrow=TRUE))
        colnames(by.item) <- c("raw_alpha","std.alpha","G6(smc)","average_r")
        rownames(by.item) <- colnames(x)
       
        Vt <- sum(R)
        item.r <- colSums(R)/sqrt(Vt)
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
        	alpha.total <- data.frame(alpha.total,mean=mean.t,sd=sdev)
        	colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)","average_r","mean","sd")
        	rownames(alpha.total) <- ""
        	stats <- data.frame(n=t.valid,r =item.r,r.cor = item.rc,r.drop = r.drop,mean=item.means,sd=item.sd)} else {
        	alpha.total <- data.frame(alpha.total)
        	        colnames(alpha.total) <- c("raw_alpha","std.alpha","G6(smc)" ,"average_r")
        	        rownames(alpha.total) <- ""

                  stats <- data.frame(r =item.r,r.cor = item.rc)
        	}
       	rownames(stats) <- colnames(x)
       	
       	
        result <- list(total=alpha.total,alpha.drop=by.item,item.stats=stats,response.freq=response.freq,keys=keys,scores = total,nvar=nvar,call=cl,title=title)
        class(result) <- c("psych","alpha")
        return(result)
     
    }
  #modified Sept 8, 2010 to add r.drop feature  
  #modified October 12, 2011 to add apply to the sd function
  #modified November 2, 2010 to use sd instead of SD
  #January 30, 2011  - added the max category parameter (max)
  #June 20, 2011 -- revised to add the check.keys option as suggested by Jeremy Miles
  #Oct 3, 2013 check for variables with no variance and drop them with a warning