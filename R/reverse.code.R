"reverse.code" <- 
function(keys,items,mini=NULL,maxi=NULL) {
 nvar <- dim(items)[2]
 items <- as.matrix(items)
 if(is.null(maxi)) {colMax <- apply(items,2,max,na.rm=TRUE)} else {colMax <- maxi}
 if(is.null(mini)) {colMin <- apply(items,2,min,na.rm=TRUE) } else {colMin <-mini}
 colAdj <- colMax+colMin
 keys.d <- diag(keys,nvar,nvar)
 items[is.na(items)] <- -9999  #a way of using matrix operations even for missing data
 reversed <- items %*% keys.d 
 keys[keys>0] <- 0
 keys <- abs(keys*colAdj)   #now we need to add the adjustments to the ones we flipped
 new <- t(keys + t(reversed))
 new[abs(new) > 999] <- NA     
 return(new) }
 
 
 
 
 
 

