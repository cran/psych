"reverse.code" <- 
function(keys,items,mini=NULL,maxi=NULL) {
 if(is.vector(items)) {nvar <- 1} else {nvar <- dim(items)[2]}
 items <- as.matrix(items)
 if(is.null(maxi)) {colMax <- apply(items,2,max,na.rm=TRUE)} else {colMax <- maxi}
 if(is.null(mini)) {colMin <- apply(items,2,min,na.rm=TRUE) } else {colMin <-mini}
 colAdj <- colMax+colMin
 
 if(length(keys) < nvar) {     
 temp <- keys 
 if(is.character(temp)) temp <- match(temp,colnames(items))
 keys <- rep(1,nvar)
 keys[temp] <- -1
 }
if(is.list(keys) | is.character(keys)){ keys <- make.keys(items,keys)
   keys <- diag(keys)}
 keys.d <- diag(keys,nvar,nvar)
 items[is.na(items)] <- -9999  #a way of using matrix operations even for missing data
 reversed <- items %*% keys.d 
 
 adj <- abs(keys*colAdj)   #now we need to add the adjustments to the ones we flipped
 adj[keys > 0] <- 0  #added December 26 
 new <- t(adj + t(reversed))
 new[abs(new) > 999] <- NA 
 colnames(new) <- colnames(items)
 colnames(new)[keys < 0] <- paste(colnames(new)[keys < 0],"-",sep="")    
 return(new) }
 
 #corrected April 3, 2010 to properly do matrix addition
#modified Sept 14, 2013  to allow symbolic names and to allow for just specifying the ones to reversed
#suggested by David Stanley
#fixed bug reported by Jian Jin (26/12/13)

 
 "rescale" <-  function(x,mean=100,sd=15,df=TRUE) {if(df) {x <- data.frame(t(t(scale(x))*sd+mean))
} else {x <- t( t(scale(x))*sd +mean)}
return(x)
}
 
 
"scrub" <- 
 function (x, where, min, max, isvalue,newvalue,cuts=NULL) 
{  
    if (missing(min))  min <- -Inf
    if (missing(max))  max <- Inf
    if (missing(isvalue))   isvalue <- Inf
    if (missing(where))  where <- 1:NCOL(x)
    maxlength <- max(length(isvalue),length(min),length(max),length(where))
    if(missing(newvalue)) newvalue <- rep(NA,maxlength)
    if (length(min) == 1)   min <- rep(min, ncol(x))
    if (length(max) == 1)  max <- rep(max, ncol(x))
    if(length(where) == 1) where <- rep(where,maxlength)
    if(length(isvalue) ==1) isvalue <- rep(isvalue,maxlength)
    if(length(newvalue) ==1) newvalue <- rep(newvalue,maxlength)
    
    if(is.null(cuts)) {for(k in 1: maxlength) {
    i <- where[k]
     x[!is.finite(x[,i]),i] <- NA
      if(is.numeric(x[,i])) {  x[(!is.na(x[, i]) & (x[, i] < min[k])), i] <- newvalue[k]
        x[(!is.na(x[, i]) & (x[, i] > max[k])), i] <- newvalue[k] 
        }
        x[(!is.na(x[, i]) & (x[, i] == isvalue[k])), i] <- newvalue[k]
       }} else {
          n.cuts <- length(cuts)
          for (j in 1:n.cuts) {
          for(k in 1:length(where)) {
           x[ (!is.na(x[where[k]])) & (x[where[k]]  >= cuts[j] ) & (x[where[k]] <  cuts[j+1]),where[k]] <- j }
                                     }
                  
     }
   
    return(x)
}
#added Sept 11, 2010
#modified December 6, 2010 to allow recoding
#modified December 3, 2011 to be more general
#modifed January 8, 2012 to be a  bit more flexible
#modified April 11, 2019 to handle character data
#fixed August 9, 2019 to correctly process is.numeric
#modified December 12th ,2020 to include the "bucketing" option
#modified December 2nd, 2021 to make Nan and inf values NA





