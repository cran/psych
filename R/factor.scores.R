"factor.scores" <- function(x,f) {
     if(!is.matrix(f)) f <- loadings(f)
     r <- cor(x,use="pairwise")   #find the correlation matrix from the data
     w <- try(solve(r,f),silent=TRUE )  #these are the factor weights
     if(class(w)=="try-error") {message("In factor.scores, the correlation matrix is singular, an approximation is used")
     ev <- eigen(r)
     ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
       r <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
       diag(r)  <- 1
      w <- solve(r,f)}
     scores <- scale(x) %*% w    #standardize the data before doing the regression
     return(scores) }
     #how to treat missing data?  see score.item
     
   