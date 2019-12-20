"Pinv" <- function(X,tol = sqrt(.Machine$double.eps)) {
   svdX <- svd(X)
   p <- svdX$d > max(tol * svdX$d[1],0 )
   if(all(p)){ Pinv <- svdX$v %*% (1/svdX$d * t(svdX$u)) } else {
    Pinv <- svdX$v[,p,drop=FALSE] %*% (1/svdX$d[p] * t(svdX$u[,p,drop=FALSE])) } 
    return(Pinv)
}