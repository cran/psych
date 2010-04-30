"describe" <-
function (x,na.rm=TRUE,interp=FALSE,skew=TRUE,ranges=TRUE,trim=.1)   #basic stats after dropping non-numeric data
                                                 #slightly faster if we don't do skews
{                      
 cl <- match.call()
#first, define a local function
    valid <- function(x) {sum(!is.na(x))}
    if(!na.rm) x <- na.omit(x)   #just complete cases
   		 
    if (is.null(dim(x)[2]))  {        #do it for vectors or 
    	    len  <- 1
    	    stats = matrix(rep(NA,9),ncol=9)    #create a temporary array
			stats[1, 1] <-  valid(x )			
			stats[1, 2] <-  mean(x, na.rm=na.rm )
			if(interp) {stats[1, 3] <- interp.median(x,na.rm=na.rm  ) }  else {stats[1,3] <- median(x,na.rm=na.rm) }
			stats[1,9] <- mean(x,na.rm=na.rm, trim=trim)
			stats[1, 4] <-  min(x, na.rm=na.rm )
			stats[1, 5] <-  max(x, na.rm=na.rm )
			stats[1, 6] <-  skew(x,na.rm=na.rm  )
			stats[1,7] <-  mad(x,na.rm=na.rm) 
			stats[1,8] <-  kurtosi(x,na.rm=na.rm) 
			vars <- 1
    	}   else  {
    	len = dim(x)[2]     #do it for matrices or data.frames 
    	
   stats = matrix(rep(NA,len*9),ncol=9)    #create a temporary array
   rownames(stats) <- colnames(x)
   stats[,1] <- apply(x,2,valid)
   vars <- c(1:len)
   ##adapted from the pairs function to convert logical or categorical to numeric 
    for(i in seq_along(x)) {
            if(is.factor(x[[i]]) || is.logical(x[[i]]))
              { x[[i]] <- as.numeric(x[[i]])
              rownames(stats)[i] <- paste(rownames(stats)[i],"*",sep="")}
            if(!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
        
   if(is.matrix(x))  {stats[,2] <- colMeans(x, na.rm=na.rm )} else {stats[,2] <- mean(x,na.rm=na.rm)}
    if (skew) {stats[, 6] <-  skew(x,na.rm=na.rm  )
               stats[,8] <- kurtosi(x,na.rm=na.rm)}
    
    for (i in 1:len) {
    	if (is.numeric(x[,i])) {   #just do this for numeric data
			
		if (ranges) {
			
			if(interp) {stats[i, 3] <- interp.median(x[,i],na.rm=na.rm  ) }  else {stats[i,3] <- median(x[,i],na.rm=na.rm) }
			stats[i,7] <-   mad(x[,i], na.rm=na.rm)
			stats[i, 4] <-  min(x[,i], na.rm=na.rm )
			stats[i, 5] <-  max(x[,i], na.rm=na.rm )
			stats[i,9]  <-  mean(x[,i],na.rm=na.rm,trim=trim)
					} #ranges
		   		 }#is.numeric
        	}# i loop	
    	} #else loop
    if (ranges)
    	{if(skew){temp <-  data.frame(n = stats[,1],mean=stats[,2], sd = SD(x), median = stats[, 
        3],trimmed =stats[,9], mad = stats[,7], min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],skew = stats[, 6], kurtosis = stats[,8])}
         
      	 else {temp <-  data.frame(n = stats[,1],mean=stats[,2], sd = SD(x), median = stats[, 
        3],trimmed =stats[,9],mad = stats[,7],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4])}}
        
        else {if(skew){temp <-  data.frame(n = stats[,1],mean=stats[,2], sd = SD(x),skew = stats[, 6], kurtosis = stats[,8])}
       else {temp <-  data.frame(n = stats[,1],mean=stats[,2], sd = SD(x))}}
                
    answer <-data.frame(var=vars,temp, se = temp$sd/sqrt(temp$n))
   
    class(answer) <- c("psych","describe","data.frame")
    return(answer)
}

