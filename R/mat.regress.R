"mat.regress" <-
function(m,x,y,digits=2)  {
 #a function to extract subsets of variables (a and b) from a correlation matrix m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
  #if there are NAs in the matrix, these are set to -999 before the reordering and then set back to NA 
  m[is.na(m)] <- -9999999
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
     	a.matrix <- reorder[1:numx,1:numx]
     	b.matrix <- reorder[1:numx,(numx+1):nxy]
     	model.mat <- solve(a.matrix,b.matrix)       #solve the equation bY~aX
     	if (length(y) >1 ) { rownames(model.mat) <- rownames(m)[x]
     	 colnames(model.mat) <- colnames(m)[y]
     	 
     	R2 <- colSums(model.mat * b.matrix) }
     	 else { R2 <- sum(model.mat * b.matrix)
     	 
     	 names(model.mat) <- rownames(m)[x]
     	 names(R2) <- colnames(m)[y]}
     	 R <- sqrt(R2)
     	mat.regress <- list(beta=round(model.mat,2),R = round(R,2),R2=round(R2,2))
     	return(mat.regress)
     	}
#modified July 12 to allow for NA in the overall matrix
