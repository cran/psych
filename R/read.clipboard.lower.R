#adapted from John Fox's read.moments function
"read.clipboard.lower" <-
function( diag = TRUE,names=NULL,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {
           xij <- scan(file("clipboard")) }
    else { xij <- scan(pipe("pbpaste"))}
    m <- length(xij)
    d <- if (diag)    1  else -1
    n <- floor((sqrt(1 + 8 * m) - d)/2) #solve the quadratic for n
    if (m != n * (n + d)/2)   stop("wrong number of elements (cannot make square matrix)")
    X <- diag(n)
    X[upper.tri(X, diag = diag)] <- xij
  diagonal <- diag(X) 
   X <- t(X) + X 
    diag(X) <- diagonal
   if(is.null(names)) names <- paste("V",1:n,sep="")
    rownames(X) <- colnames(X) <- names
    return(X)
   }


"read.clipboard.upper" <-
function( diag = TRUE,names=NULL,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {
           xij <- scan(file("clipboard")) }
    else { xij <- scan(pipe("pbpaste"))}
    m <- length(xij)
    d <- if (diag)    1  else -1
    n <- floor((sqrt(1 + 8 * m) - d)/2) #solve the quadratic for n
    if (m != n * (n + d)/2)   stop("wrong number of elements (cannot make square matrix)")
    X <- diag(n)
    X[lower.tri(X, diag = diag)] <- xij
  diagonal <- diag(X) 
   X <- t(X) + X 
    diag(X) <- diagonal
   if(is.null(names)) names <- paste("V",1:n,sep="")
    rownames(X) <- colnames(X) <- names
    return(X)
   }
   
   
  