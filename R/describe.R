"describe" <-
function (x, digits = 2,na.rm=TRUE,skew=FALSE,ranges=TRUE)   #basic stats after dropping non-numeric data
                                                 #much faster if we don't do skews
{                         #first, define a local function
    valid <- function(x) {sum(!is.na(x))}
   		 
    if (is.vector(x) )  {        #do it for vectors or 
    	    len  <- 1
    	    stats = matrix(rep(NA,7),ncol=7)    #create a temporary array
			stats[1, 1] <-  valid(x )			
			stats[1, 2] <-  mean(x, na.rm=na.rm )
			stats[1, 3] <-  median(x,na.rm=na.rm  )
			stats[1, 4] <-  min(x, na.rm=na.rm )
			stats[1, 5] <-  max(x, na.rm=na.rm )
			stats[1, 6] <-  skew(x,na.rm=na.rm  )
			stats[1,7] <-  mad(x,na.rm=na.rm) 
			
    	}   else  {
    	len = dim(x)[2]     #do it for matrices or data.frames 
    	
   stats = matrix(rep(NA,len*7),ncol=7)    #create a temporary array
   stats[,1] <- apply(x,2,valid)
   stats[,2] <- colMeans(x, na.rm=na.rm )
    if (skew) {stats[, 6] <-  skew(x,na.rm=na.rm  )}
    
    for (i in 1:len) {
    	if (is.numeric(x[,i])) {   #just do this for numeric data
			
		if (ranges) {
			stats[i, 3] <-  median(x[,i],na.rm=na.rm  )	
			stats[i,7] <- mad(x[,i], na.rm=na.rm)
			stats[i, 4] <-  min(x[,i], na.rm=na.rm )
			stats[i, 5] <-  max(x[,i], na.rm=na.rm )
					} #ranges
		   		 }#is.numeric
        	}# i loop	
    	} #else loop
    if (ranges)
    	{if(skew){temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE), median = stats[, 
        3],mad = stats[,7], min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4],skew = stats[, 6])}
         
      	 else {temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE), median = stats[, 
        3],mad = stats[,7],min= stats[,4],max=stats[,5], range=stats[,5]-stats[,4])}}
        
        else {if(skew){temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE),skew = stats[, 6])}
       else {temp <-  data.frame(var = seq(1:len),n = stats[,1],mean=stats[,2], sd = sd(x,na.rm=TRUE))}}
                
    answer <-  round(data.frame(temp, se = temp$sd/sqrt(temp$n)),  digits)
     return(answer)
}

