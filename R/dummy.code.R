"dummy.code" <- function(x) {
 t <- table(x)
 lt <- length(t)
 n.obs <- length(x)
 new <- matrix(0,nrow=n.obs,ncol=lt)
 xlev <- as.factor(x)
 for (i in 1:n.obs) {new[i,xlev[i]] <- 1}
 colnames(new) <- names(t)
 return(new)
}