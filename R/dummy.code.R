"dummy.code" <- function(x,group=NULL,na.rm=TRUE) {
t <- table(x)
 lt <- length(t)
 n.obs <- length(x)
 if(is.null(group)) {new <- matrix(0,nrow=n.obs,ncol=lt) #the case of no grouping information
 if(na.rm) {new[is.na(x),] <- NA } #added 10/20/17
 xlev <- as.factor(x)
 for (i in 1:n.obs) {
new[i,xlev[i]] <- 1} 
 
 colnames(new) <- names(t) 
 }  else {new <- rep(0,n.obs)  #the alternative is to combine categories 
  xlev <- as.factor(x)
   if(na.rm) {new[is.na(x)] <- NA } 
 for (i in 1:n.obs) {
new[i] <- xlev[i] %in% group} 
 }
 
 return(new)
}
