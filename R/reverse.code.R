"reverse.code" <- 
function(keys,items,mini=NULL,maxi=NULL) {
 if(is.vector(items)) {nvar <- 1} else {nvar <- dim(items)[2]}
 items <- as.matrix(items)
 if(is.null(maxi)) {colMax <- apply(items,2,max,na.rm=TRUE)} else {colMax <- maxi}
 if(is.null(mini)) {colMin <- apply(items,2,min,na.rm=TRUE) } else {colMin <-mini}
 colAdj <- colMax+colMin
 keys.d <- diag(keys,nvar,nvar)
 items[is.na(items)] <- -9999  #a way of using matrix operations even for missing data
 reversed <- items %*% keys.d 
 keys[keys > 0] <- 0
 keys <- abs(keys*colAdj)   #now we need to add the adjustments to the ones we flipped
 new <- t(keys + t(reversed))
 new[abs(new) > 999] <- NA     
 return(new) }
 
 
 "rescale" <-  function(x,mean=100,sd=15,df=TRUE) {if(df) {x <- data.frame(t(t(scale(x))*sd+mean))
} else {x <- t( t(scale(x))*sd +mean)}
return(x)
}
#corrected April 3, 2010 to properly do matrix addition
 
 
 
 
"scrub" <- 
function(x,where,min,max,isvalue) {
if(missing(min)) min <- -Inf
if(missing(max)) max <-  Inf
if(missing(isvalue)) isvalue <- Inf
if(missing(where)) where <- 1:dim(x)[2]
if(length(min) == 1) min <- rep(min,(length(where)))	
if(length(max) == 1) max <- rep(max,(length(where)))	
if(length(isvalue) == 1) isvalue <- rep(isvalue,(length(where)))	
k <- 1
for(i in where) {
x[(!is.na(x[,i]) &(x[,i] < min[k])),i]<- NA
x[(!is.na(x[,i]) &(x[,i] > max[k])),i]<- NA
x[(!is.na(x[,i]) &(x[,i] == isvalue[k])),i]<- NA
k <- k+1}
return(x)}
#added Sept 11, 2010

#does not work
"scrubber" <- 
function(x,who,where,min,max,isvalue) {
if(missing(who)) who <- 1:dim(x)[1]
if(missing(min)) min <- -Inf
if(missing(max)) max <-  Inf
if(missing(isvalue)) isvalue <- Inf
if(missing(where)) where <- 1:dim(x)[2]
if(length(min) == 1) min <- rep(min,(length(where)))	
if(length(max) == 1) max <- rep(max,(length(where)))	
if(length(isvalue) == 1) isvalue <- rep(isvalue,(length(where)))	
k <- 1
for(i in where) {
x[(!is.na(x[who,i]) &&(who[x[who,i] < min[k]])),i]<- NA
x[(!is.na(x[who,i]) &&(x[who,i] > max[k])),i]<- NA
#x[(!is.na(x[who,i]) &&(x[who,i] == isvalue[k])),i]<- NA
k <- k+1}
return(x)}
#added Sept 11, 2010

