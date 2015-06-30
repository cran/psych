
#added 1/11/14
#modified 12/29/14 to handle cases with non-numeric data
#for fast descriptions
describeData <- 
function (x, head = 4, tail = 4) 
{
    valid <- function(x) {
        sum(!is.na(x))
    }
    nvar <- ncol(x)
    all.numeric <- nvar
    ans <- matrix(NA,nrow=nvar,ncol=2)
    nobs <- nrow(x)
    cc <- 0
    cc <- try(complete.cases(x),silent=TRUE) 
    if(class(cc) == "try-error") cc <- NA
    cc <- sum(cc,na.rm=TRUE)
   for (i in 1:nvar) {
        
        if (is.numeric(x[,i])) {ans[i,2] <- 1 } else {
            if ((is.factor(x[,i])) || (is.logical(x[,i]))) {
               ans[i,2]  <- 2
            } else {
                if (is.character(x[,i])) {
               ans[i,2] <- 3
                } else {ans[i,2] <- 4}
            }
        }
        ans[i,1]  <- valid(x[,i])
    }

  
    if (is.numeric(unlist(x))) {
        all.numeric <- TRUE
    }
    else {
        all.numeric <- FALSE
    }
    H1 <- t(x[1:head,1:nvar])
    T1 <- t(x[(nobs-tail+1):nobs,1:nvar])
    temp <- data.frame(V=1:nvar,ans,H1,T1)
    
    colnames(temp) <- c("variable #", "n.obs", "type", paste("H", 
        1:head, sep = ""), paste("T", 1:tail, sep = ""))
   rownames(temp)[temp[,"type"]!=1] <- paste(rownames(temp)[temp[,"type"]!=1],"*",sep="")
    result <- (list(n.obs = nobs, nvar = nvar, all.numeric = all.numeric, 
        complete.cases = cc, variables = temp))
    class(result) <- c("psych", "describeData")
    return(result)
}


#changed October 12, 2011 to use apply because mean and sd are deprecated for data.frames
#modified January 10, 2014 to add the check option to improve speed.  A few other improvements
#modified December 2014 to add the fast option for large data sets
#modified May 21, 2015 to allow non-numeric data to be described (but with a warning)
"describe" <-
function (x,na.rm=TRUE,interp=FALSE,skew=TRUE,ranges=TRUE,trim=.1,type=3,check=TRUE,fast=NULL)   #basic stats after dropping non-numeric data
                                                 #slightly faster if we don't do skews
{                      
 cl <- match.call()
#first, define a local function
    valid <- function(x) {sum(!is.na(x))}
    if(!na.rm) x <- na.omit(x)   #just complete cases
   	
   	if(is.null(fast)) {
   	    if (prod(dim(x)) > 10^7) {fast <- TRUE } else {fast <- FALSE}}  #the default is to use fast for large data sets
   	if(fast) {skew <- FALSE
   	         }
   		 
    if ( is.null(dim(x)[2]))  {        #do it for vectors or 
    	    len  <- 1
    	    stats = matrix(rep(NA,10),ncol=10)    #create a temporary array
			stats[1, 1] <-  valid(x )			
			stats[1, 2] <-  mean(x, na.rm=na.rm )
			stats[1,10] <- sd(x,na.rm=na.rm)
			if(interp) {stats[1, 3] <- interp.median(x,na.rm=na.rm  ) }  else {stats[1,3] <- median(x,na.rm=na.rm) }
			stats[1,9] <- mean(x,na.rm=na.rm, trim=trim)
			stats[1, 4] <-  min(x, na.rm=na.rm )
			stats[1, 5] <-  max(x, na.rm=na.rm )
			stats[1, 6] <-  skew(x,na.rm=na.rm,type=type  )
			stats[1,7] <-  mad(x,na.rm=na.rm) 
			stats[1,8] <-  kurtosi(x,na.rm=na.rm,type=type) 
			vars <- 1
    	} else {
    	
   stats = matrix(rep(NA,ncol(x)*10),ncol=10)    #create a temporary array
   rownames(stats) <- colnames(x)
   stats[,1] <- apply(x,2,valid)
   vars <- c(1:ncol(x))
   ##adapted from the pairs function to convert logical or categorical to numeric 
    if(check) {
    for(i in 1:ncol(x)) {   
        if(!is.numeric(x[[i]] ))  {
                                 if(fast)  {x[[i]] <- NA} else {
                                  if(is.factor(unlist(x[[i]]))) {   #fixed 5/21/15
                                              x[[i]] <- as.numeric(x[[i]]) 
                          } else {x[[i]] <- NA} }
               
              rownames(stats)[i] <- paste(rownames(stats)[i],"*",sep="")}
              } }
      
    x <- as.matrix(x)
   
    stats[,2] <- apply(x, 2,mean,na.rm=na.rm )
    stats[,10] <- apply(x,2,sd,na.rm=na.rm)
    
    
    if (skew) {stats[, 6] <-  skew(x,na.rm=na.rm,type=type  )
               stats[,8] <- kurtosi(x,na.rm=na.rm,type=type)}
    
    
     if(ranges) {
            if(fast) {
          stats[,4] <- apply(x,2,min,na.rm=na.rm)
          stats[,5] <- apply(x,2,max,na.rm = na.rm)
          } else {
     			stats[, 4] <-  apply(x,2,min, na.rm=na.rm )
			    stats[, 5] <-  apply(x,2,max, na.rm=na.rm )
    		    stats[,7] <-   apply(x,2,mad, na.rm=na.rm)
			    stats[,9]  <- apply(x,2, mean,na.rm=na.rm,trim=trim)
		if(interp) {stats[, 3] <- apply(x,2,interp.median,na.rm=na.rm  ) }  else {stats[,3] <- apply(x,2,median,na.rm=na.rm) } 
		}} 
		
    
    
    }  #end of maxtrix input
     #now summarize the results 
    if (ranges) {
    			if(skew){
    			   answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 3],trimmed =stats[,9], mad = stats[,7], min= stats[,4],max=stats[,5],
    			   range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1]))  #the typical (maximum) case
                        } else {
                                 if(!fast) {
                                     answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[,3],trimmed =stats[,9],
                                     mad = stats[,7],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/sqrt(stats[,1]))  #somewhat shorter 
                                           } else {answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/sqrt(stats[,1]))  #even shorter 
                                                   }
                               }} else {
                                       if(skew){answer <-  data.frame(vars,n = stats[,1],mean=stats[,2], sd =stats[,10],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1]))} else {
               answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],se=stats[,10]/sqrt(stats[,1]))}}  #the minimal case -- fixed 1/2/2014
                
   # answer <-data.frame(var=vars,temp, se = temp$sd/sqrt(temp$n))  #replaced with forming the answer in the if statements  10/1/2014 to improve memory management
   
    class(answer) <- c("psych","describe","data.frame")
    return(answer)
}
