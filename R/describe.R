	#Added 2/10/24
	#This is a parellel version of describe  --much faster
	#just calls the old describe and parallelizes it
	
	describe <- function(x,na.rm=TRUE,interp=FALSE,skew=TRUE,ranges=TRUE,trim=.1,type=3,check=TRUE,fast=NULL,quant=NULL,
    IQR=FALSE,omit=FALSE,data=NULL,size=50)  {
    if(inherits(x,"formula")) {ps <- fparse(x)   #group was specified, call describeBy
	if(missing(data)) { 
	 x <- get(ps$y)
	group <- x[,ps$x]} else {x <- data[ps$y]
	group <- data[ps$x]
	}
	describeBy(x,group=group,na.rm=na.rm,interp=interp,skew=skew,ranges=ranges,trim=trim,type=type,check=check,fast=fast,quant=quant,IQR=IQR,omit=omit,data=data)
  } else {                   
 cl <- match.call()
#first, define a local function
    valid <- function(x) {sum(!is.na(x))}
    if(!na.rm) x <- na.omit(x)   #just complete cases
   	 if(is.null(fast)) {   #don't reset fast inside describe.1 loop
   	   if (prod(dim(x)) > 10^7) {fast <- TRUE } else {fast <- FALSE}
   	   } 
nvar <- NCOL(x)
v.names <- colnames(x)
if(nvar <  size) {result <- describe.1(x=x,na.rm=na.rm, interp=interp, skew=skew, ranges=ranges,
                             trim=trim, type=type,check=check,fast=fast,quant=quant,
                             IQR = IQR, omit=omit, data=data)} else {
                             
    n.steps <- ceiling(nvar/size)
   
   short <- function(i) {
   	loweri <- (i-1) * size + 1
   	upperi <- min(i * size,nvar)
     res <- describe.1(x=x[,loweri:upperi],na.rm=na.rm, interp=interp, skew=skew, ranges=ranges,
                             trim=trim, type=type,check=check,fast=fast,quant=quant,
                             IQR = IQR, omit=omit, data=data)   
                   return(res)                 
                             } 

        result <- mcmapply(short,c(1:n.steps))
        n.result <- NROW(result)   
                                   
       names.result <- rownames(result)
         
        temp <- NULL
       for( i in 1:NCOL(result)) {tt <- matrix(unlist(result[,i]),ncol=n.result)  #fix to be dynamic
           # tt <- result[[i]]
           temp <- rbind(temp,tt)}
        
        # rownames(temp) <- colnames(x) 
        # colnames(temp) <-  c("vars" ,   "n" ,  "mean",     "sd", "median",  "min" ,   "max", "range" ,    "se") #fix to be dynamic
        colnames(temp) <- names.result
         temp[,1 ] <-1:nvar
         rownames(temp )  <- v.names
      result <- as.data.frame(temp)}
      
      class(result) <- c("psych","Pdescribe","data.frame")
      return(result)
} 
}



"describeFast" <- function(x) {
 if(!inherits(x[1], "data.frame")) x <- fix.dplyr(x)    #to get around a problem created by dplyr
nvar <- NCOL(x)
nobs <- NROW(x)
valid <- colSums(!is.na(x))
temp <- matrix(NA,nrow=nvar,ncol=4)
for(i in 1:nvar) {temp[i,1] <- is.numeric(x[1,i])
                temp[i,2] <- is.factor(x[1,i])
                temp[i,3] <- is.logical(x[1,i])
                temp[i,4] <- is.character(x[1,i]) }
ttt <- which(temp[,1:4] == TRUE,arr.ind=TRUE)
if(nvar > 1) {
ttt <- dfOrder(ttt,"row")}
temp <- cbind(temp,ttt["col"])
colnames(temp) <- c("numeric","factor","logical","character","type")

 cc <- try(complete.cases(x),silent=TRUE) 
    if(inherits(cc,"try-error"))  cc <- NA
    cc <- sum(cc,na.rm=TRUE)
all.numeric <- sum(temp[,1])
all.factor <- sum(temp[,2])

result.df <- data.frame(var=1:nvar,n.obs=valid,temp)
result<- list(nvar=nvar,n.obs =nobs,complete.cases = cc,numeric=all.numeric,factors=all.factor,result.df=result.df)
class(result) <- c("psych","describeFast")
return(result)
}

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
    if(inherits(cc, "try-error")) cc <- NA
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
#further modified June 21, 2016 to allow for character input as well as well reporting quantiles 
#tried to improve the speed by using multicores, but this requires using s or lapply which don't do what I need.
#Added the as.factor  to the convert to numeric for non-numeric data
#Improved July 9, 2020 to allow formula input for grouping variables
#Further improved Sept 02, 2020 to allow formula input with data specified to match other formula features
#include median with fast output (has a very small slowing effect)
"describe.1" <-
function (x,na.rm=TRUE,interp=FALSE,skew=TRUE,ranges=TRUE,trim=.1,type=3,check=TRUE,fast=NULL,
        quant=NULL,IQR=FALSE,omit=FALSE,data=NULL)   #basic stats after dropping non-numeric data
                             #slightly faster if we don't do skews
{ 
if(inherits(x,"formula")) {ps <- fparse(x)   #group was specified, call describeBy
	if(missing(data)) { 
	 x <- get(ps$y)
	group <- x[,ps$x]} else {x <- data[ps$y]
	group <- data[ps$x]
	}
	describeBy(x,group=group,na.rm=na.rm,interp=interp,skew=skew,ranges=ranges,trim=trim,type=type,check=check,fast=fast,quant=quant,IQR=IQR,omit=omit,data=data)
  } else {                   
 cl <- match.call()
#first, define a local function
    valid <- function(x) {sum(!is.na(x))}
    if(!na.rm) x <- na.omit(x)   #just complete cases
   	
   	#if(is.null(fast)) {   #don't reset fast inside describe.1 loop
   	 #   if (prod(dim(x)) > 10^7) {fast <- TRUE } else {fast <- FALSE}}  #the default is to use fast for large data sets
   #	if(fast) {skew <- FALSE
   	#         }
   	numstats <- 10 + length(quant)	+ IQR 
   	cn <- "X1"
    if ( NCOL(x) < 2)  {if(is.data.frame(x)) {      #
                if( !is.numeric(x[,1])) {warning ("You were trying to describe a non-numeric data.frame or vector which describe converted  to numeric.")
                    x[,1] <- as.numeric(as.factor(x[,]))
                    } 
              if(!is.null(colnames(x))) {cn <- colnames(x)} else {cn  <- "X1"}
             
               x <- x[,1] }   #getting around the problem of single column data frames       
                #do it for vectors or 
    	    len  <- 1
    	    nvar <- 1
    	    stats = matrix(rep(NA,numstats),ncol=numstats)    #create a temporary array
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
			 if(!is.null(quant)) { Qnt <- quantile(x,prob=quant,na.rm=TRUE)
             stats[1,(IQR+11):numstats] <- t(Qnt)}
             if(IQR) {Quart <- t(quantile(x,prob=c(.25,.75),na.rm=TRUE))
             Iqr <- Quart[,2] -Quart[,1] 
             stats[1,11] <- Iqr
             }
            
			rownames(stats) <- cn
    	} else {
   nvar <- ncol(x)	
   stats = matrix(rep(NA,nvar*numstats),ncol=numstats)    #create a temporary array
   
   if(is.null(colnames(x))) colnames(x) <- paste0("X",1:ncol(x))
   rownames(stats) <- colnames(x)
   stats[,1] <- apply(x,2,valid)
   vars <- c(1:nvar)
   ##adapted from the pairs function to convert logical or categorical to numeric 
   select <- 1:nvar

    if(!is.matrix(x) && check) {  #does not work for matrices
    for(i in 1:nvar) {   
        if(!is.numeric(x[[i]] ))  {
                                 if(fast)  {x[[i]] <- NA} else {
                                 if(omit) {select[i] <- NA}
                                  if(is.factor(unlist(x[[i]])) | is.character(unlist(x[[i]]))) {  x[[i]] <- as.numeric(as.factor(x[[i]])) 
                                   rownames(stats)[i] <- paste(rownames(stats)[i],"*",sep="")
                                  
                          } else {x[[i]] <- NA} 
               
             }
              }
              } 
             } 
             
              
    select <- select[!is.na(select)]
    
    x <- as.matrix(x[,select])
    vars <- vars[select]
   
    stats <- stats[select,]
    if(!is.numeric(x)) {message("Converted non-numeric matrix input to numeric.  Are you sure you wanted to do this. Please check your data")
                x <- matrix(as.numeric(x),ncol=nvar)
               rownames(stats) <- paste0(rownames(stats),"*")} 
   
    stats[,2] <- apply(x, 2,mean,na.rm=na.rm )
    stats[,10] <- apply(x,2,sd,na.rm=na.rm)
    
    
    if (skew) {stats[, 6] <-  skew(x,na.rm=na.rm,type=type  )
               stats[,8] <- kurtosi(x,na.rm=na.rm,type=type)}
    
    
     if(ranges) {
            if(fast) {
          stats[,4] <- apply(x,2,min,na.rm=na.rm)
          stats[,5] <- apply(x,2,max,na.rm = na.rm)
          if(interp) {stats[, 3] <- apply(x,2,interp.median,na.rm=na.rm  ) }  else {stats[,3] <- apply(x,2,median,na.rm=na.rm) } 
          } else {
     			stats[, 4] <-  apply(x,2,min, na.rm=na.rm )
			    stats[, 5] <-  apply(x,2,max, na.rm=na.rm )
    		    stats[,7] <-   apply(x,2,mad, na.rm=na.rm)
			    stats[,9]  <- apply(x,2, mean,na.rm=na.rm,trim=trim)
			if(interp) {stats[, 3] <- apply(x,2,interp.median,na.rm=na.rm  ) }  else {stats[,3] <- apply(x,2,median,na.rm=na.rm) } 
		}} 
		
    if(!is.null(quant)) { Qnt <- apply(x,2,quantile,prob=quant,na.rm=TRUE)
    stats[,(IQR+11):numstats] <- t(Qnt)}
    
    if(IQR) {Quart <- t(apply(x,2,quantile,prob=c(.25,.75),na.rm=TRUE))
             Iqr <- Quart[,2] - Quart[,1] 
             stats[,11] <- Iqr
             }
    }  #end of maxtrix input
     #now summarize the results 
     if (numstats > (10 + IQR)) {
     colnames(stats)[(11+IQR):numstats] <- paste0("Q",quant[1:length(quant)])}
    
    #the following output was cleaned up on June 22, 2016 added the quantile information.
    
    #the various options are ranges, skew, fast, numstats > 10
    if(fast) { answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],median=stats[,3],se=stats[,10]/sqrt(stats[,1])) }  #minimal case
    
    #if((!skew) && ranges) {answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/sqrt(stats[,1])) }
    if(skew) {
        if(ranges) {if(fast) { answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 3], min= stats[,4],max=stats[,5],
    				 range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1])) } else {
    				  answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 3],trimmed =stats[,9], mad = stats[,7], min= stats[,4],max=stats[,5],
    				 range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1])) }
    			   } else {
            		 answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],skew = stats[, 6], 
            		        kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1])) }
             } else {if(ranges) {answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],median=stats[,3],
             min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/sqrt(stats[,1])) } else {
             answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],se=stats[,10]/sqrt(stats[,1]))   }
            }
    if(IQR) answer <- data.frame(answer,IQR=stats[,11])
          
    if (numstats > (10+ IQR)) {if(nvar > 1 ) {answer <- data.frame(answer, stats[,(IQR+11):numstats])   #add the quantile information
                                     } else { 

                                       answer <- data.frame(answer, t(stats[,(IQR+11):numstats])) 
                                       rownames(answer) <- cn}
                        }     
                         
    
#  {if (ranges) {
#     	if(skew){ 
#     		 
#     		 
#     	   if(numstats > 10) {  answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 3],trimmed =stats[,9], mad = stats[,7], min= stats[,4],max=stats[,5],
#     			   range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8], se=stats[,10]/sqrt(stats[,1]), stats[,(11:numstats)])}  else {
#     		 answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 3],trimmed =stats[,9], mad = stats[,7], min= stats[,4],max=stats[,5],
#     		 range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1]))}  #the typical (maximum) case
#                            } else {
#                                  if(!fast) {
#                                      if(numstats > 10) {
#                                      answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[,3],trimmed =stats[,9],
#                                      mad = stats[,7],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/sqrt(stats[,1]))  #somewhat shorter 
#                                            } else {answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/sqrt(stats[,1]))  #even shorter 
#                                                    }
#                                }} else {
#                                        if(skew){answer <-  data.frame(vars,n = stats[,1],mean=stats[,2], sd =stats[,10],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/sqrt(stats[,1]))} else {
#                answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],se=stats[,10]/sqrt(stats[,1]))}}  #the minimal case -- fixed 1/2/2014
#      }         
# answer <-data.frame(var=vars,temp, se = temp$sd/sqrt(temp$n))  #replaced with forming the answer in the if statements  10/1/2014 to improve memory management

    class(answer) <- c("psych","Pdescribe","data.frame")
    return(answer)
    }
}
#minor patch on 9/10/23 for case of single variables with quants