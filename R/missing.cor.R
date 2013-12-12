"missing.cor" <- function(x) {
  nvar <- ncol(x)
  Mux <- apply(x,2,mean,na.rm=TRUE)
  Varx <- apply(x,2,var,na.rm=TRUE)
  X <- scale(x,scale=FALSE) 
  Covx <- diag(Varx,ncol=nvar) 
  N <-  t(!is.na(x)) %*% (!is.na(x)) 
  
 
     
for(i in 2:nvar) {
  for (j in 1:(i-1)) {
     Covx[i,j]  <- sum(X[i]*X[j],na.rm=TRUE)
     }
     }
Covx <- Covx/(N-1)
}