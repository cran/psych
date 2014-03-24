"partial.r" <-
function(m,x,y)  {
 cl <- match.call()
   if(dim(m)[1]!=dim(m)[2]) {n.obs <- dim(m)[1]
                    m <- cor(m,use="pairwise") }
  if(!is.matrix(m)) m <- as.matrix(m)

    #first reorder the matrix to select the right variables
         nm <- dim(m)[1]
        t.mat <- matrix(0,ncol=nm,nrow=nm)
        xy <- c(x,y)
         numx <- length(x)
     	numy <- length(y)
        nxy <- numx+numy
        for (i in 1:nxy) {
     	t.mat[i,xy[i]] <- 1 }
     	
     	reorder <- t.mat %*% m %*% t(t.mat)
     	reorder[abs(reorder)  > 1] <- NA    #this allows us to use the matrix operations to reorder and pick
     	X <- reorder[1:numx,1:numx]
     	Y <- reorder[1:numx,(numx+1):nxy]
     	
        phi <- reorder[(numx+1):nxy,(numx+1):nxy] 
        phi.inv <- solve(phi)
        X.resid <- X - Y %*% phi.inv %*% t(Y)
       # sd <- diag(sqrt(1/diag(X.resid)))
       # X.resid <- sd %*% X.resid %*% sd
       X.resid <- cov2cor(X.resid)  
       
        colnames(X.resid) <- rownames(X.resid) <- colnames(m)[x]
        
        class(X.resid)  <- c("psych","partial.r")
        return(X.resid)
     	}
     #modified March 23 to use cov2cor instead of the sd line.  This makes the diagonal exactly 1.
     

