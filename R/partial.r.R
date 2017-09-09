"partial.r" <-
function(data,x,y,use="pairwise",method="pearson")  {
 cl <- match.call()
   if(!isCorrelation(data)) {n.obs <- dim(data)[1]
                    data <- cor(data,use=use,method=method) }
   m <- as.matrix(data)
 if(missing(x) & missing(y)) {X.resid <- -(solve(m))  #this is thus the image covariance matrix
    diag(X.resid) <- 1/(1- smc(m))    #adjust the diagonal to be 1/error
    X.resid <- cov2cor(X.resid)} else {
                
        xy <- c(x,y)                     
     	X <- m[x,x]
     	Y <- m[x,y]
     	phi <- m[y,y]
        phi.inv <- solve(phi)
        X.resid <- X - Y %*% phi.inv %*% t(Y)
        X.resid <- cov2cor(X.resid) 
        class(X.resid)  <- c("psych","partial.r") }
       
        return(X.resid)
     	}
     #modified March 23 to use cov2cor instead of the sd line.  This makes the diagonal exactly 1.
     #05/08/17  Completely rewritten to be easier to use and follow for the case of complete partials 

