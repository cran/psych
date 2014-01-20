#changed October 12, 2011 to use apply because mean and sd are deprecated for data.frames
#modified January 10, 2014 to add the check option to improve speed.  A few other 
"describe" <-
function (x,na.rm=TRUE,interp=FALSE,skew=TRUE,ranges=TRUE,trim=.1,type=3,check=TRUE)   #basic stats after dropping non-numeric data
                                                 #slightly faster if we don't do skews
{                      
 cl <- match.call()
#first, define a local function
    valid <- function(x) {sum(!is.na(x))}
    if(!na.rm) x <- na.omit(x)   #just complete cases
   		 
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
            #if(is.factor(x[[i]]) || is.logical(x[[i]]))    #basically this is the same question as 
            if(!is.numeric(x[[i]] ))  {
             x[[i]] <- as.numeric(x[[i]])
              rownames(stats)[i] <- paste(rownames(stats)[i],"*",sep="")}
            if(!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'describe'")
        }
        }
      
    x <- as.matrix(x)
    
    stats[,2] <- apply(x, 2,mean,na.rm=na.rm )
    stats[,10] <- apply(x,2,sd,na.rm=na.rm)
    
    if (skew) {stats[, 6] <-  skew(x,na.rm=na.rm,type=type  )
               stats[,8] <- kurtosi(x,na.rm=na.rm,type=type)}
    
    
     if(ranges) {
    		stats[,7] <-   apply(x,2,mad, na.rm=na.rm)
			stats[, 4] <-  apply(x,2,min, na.rm=na.rm )
			stats[, 5] <-  apply(x,2,max, na.rm=na.rm )
			stats[,9]  <- apply(x,2, mean,na.rm=na.rm,trim=trim)
		if(interp) {stats[, 3] <- apply(x,2,interp.median,na.rm=na.rm  ) }  else {stats[,3] <- apply(x,2,median,na.rm=na.rm) } }
    
    
    }  #end of maxtrix input
     #now summarize the results 
    if (ranges)
    	{if(skew){answer  <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 
        3],trimmed =stats[,9], mad = stats[,7], min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/stats[,1])  #the typical (maximum) case
         } else {
          answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10], median = stats[, 
        3],trimmed =stats[,9],mad = stats[,7],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],se=stats[,10]/stats[,1])}  #somewhat shorter 
        
       } else {if(skew){answer <-  data.frame(vars,n = stats[,1],mean=stats[,2], sd =stats[,10],skew = stats[, 6], kurtosis = stats[,8],se=stats[,10]/stats[,1])}
       else {answer <-  data.frame(vars=vars,n = stats[,1],mean=stats[,2], sd = stats[,10],se=stats[,10]/stats[,1])}}  #the minimal case
                
   # answer <-data.frame(var=vars,temp, se = temp$sd/sqrt(temp$n))  #replaced with forming the answer in the if statements  10/1/2014 to improve memory management
   
    class(answer) <- c("psych","describe","data.frame")
    return(answer)
}

#added 1/11/14
"describeData" <- function(x,head=4,tail=4) {
valid <- function(x) {sum(!is.na(x))}
nvar <- ncol(x)
all.numeric <- nvar
ans <- c(NA,(3+head + tail))

nobs <- nrow(x)
cc <- sum(complete.cases(x))
#rownames(ans) <- colnames(x)
temp <- apply(x,2, function(xx) {
if(!is.numeric(xx)) {
   if( (is.factor(xx)) || (is.logical(xx))) {ans[3] <- "factor/logical"} else {
   if(is.character(xx)) {ans[3] <- "character"} 
  }
}
ans[2] <- valid(xx)
ans[4:(3+head)] <- xx[1:head]
ans[(4+head):(3+head+tail)] <- xx[((nobs-tail+1):nobs)] 
ans})

temp[1,] <- 1:nvar
temp <- t(temp)
if(is.numeric(temp)) {all.numeric<- TRUE} else {all.numeric <- FALSE}
colnames(temp) <- c("variable #","n.obs", "type" ,paste("H",1:head,sep=""),paste("T",1:tail,sep=""))
result <- (list(n.obs= nobs,nvar=nvar,all.numeric=all.numeric, complete.cases = cc,variables = temp))
class(result) <- c('psych','describeData')
return(result)
}
 


