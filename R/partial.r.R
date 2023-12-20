"partial.r" <-
function(data,x,y,use="pairwise",method="pearson",part=FALSE)  {
 cl <- match.call()
 #convert formula input into prior format
 #x are the set from which y is partialled
 
if(!missing(x))  { if(inherits(x,"formula")) {
   ps <- fparse(x)
   y <- ps$y
   x <- ps$x
   z <- ps$z   #do we have any variable to partial out
   ex <- ps$ex
   
   #now put it into old form
   x <- c(y,x)
   yy <- y  #we need to keep the name(s) of the dependent variables
   y <- z   #we partial these from x 
}
}
   if(!isCorrelation(data)) {n.obs <- dim(data)[1]
    if(!missing(x) & !missing(y)) {if(!is.character(x) ) x <- colnames(data)[x]
       if(!is.character(y) ) y <- colnames(data)[y]
        data <- cor(data[,c(x,y)],use=use,method=method)
         } else {if(is.null(dim(data))) stop("Specify the rows for data (use , for all rows)")
                    data <- cor(data,use=use,method=method) }}
   m <- as.matrix(data)
 #if(missing(x) & missing(y)) {X.resid <- -(solve(m))  #this is thus the image covariance matrix
 if(missing(x) & missing(y)) {X.resid <- -(Pinv(m))  #this is thus the image covariance matrix
    diag(X.resid) <- 1/(1- smc(m))    #adjust the diagonal to be 1/error
    X.resid <- cov2cor(X.resid)
    rownames(X.resid) <-colnames(X.resid) <- colnames(m)} else {
        if(missing(x)){ x <- colnames(data)
        y <- as.character(y)
        x <- x[!x %in% y]}
               
        xy <- c(x,y)                     
     	X <- m[x,x]
     	Y <- m[x,y]
     	phi <- m[y,y]
       # phi.inv <- solve(phi)
       phi.inv <- Pinv(phi)
        X.resid <- X - Y %*% phi.inv %*% t(Y)
        
      
        if(part) if(length(yy) >1) {diag(X.resid[yy,yy]) <- 1} else {X.resid[yy,yy] <- 1}
        X.resid <- cov2cor(X.resid) 
        class(X.resid)  <- c("psych","partial.r", "matrix") }
       
        return(X.resid)
     	}
     #modified March 23 to use cov2cor instead of the sd line.  This makes the diagonal exactly 1.
     #05/08/17  Completely rewritten to be easier to use and follow for the case of complete partials 
     #modified 03/19/19 to just choose the items to correlate instead of entire matrix
     #modified 07/25/20 to use the Pseudo Inverse so that in cases of improper matrices, we still give a partial
     #modified 06/09/21 to add the matrix class to the object.
     #modified 12/3/21 to add formula input option
     #modified 12/5/23 to add the ability to do part correlations (suggested by Rick Zinbarg)
     
     
