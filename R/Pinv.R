"Pinv" <- function(R,tol = sqrt(.Machine$double.eps)) {
   svdX <- svd(R)
   p <- svdX$d > max(tol * svdX$d[1],0 )
   if(all(p)){ Pinv <- svdX$v %*% (1/svdX$d * t(svdX$u)) } else {
    Pinv <- svdX$v[,p,drop=FALSE] %*% (1/svdX$d[p] * t(svdX$u[,p,drop=FALSE])) } 
    return(Pinv)
}

"matSqrt" <- function(R,tol = sqrt(.Machine$double.eps)) {
    svdX <- svd(R)
   p <- svdX$d > max(tol * svdX$d[1],0 )
   if(all(p)){ result <- svdX$v %*% (sqrt(svdX$d) * t(svdX$u)) } else {
    result <- svdX$v[,p,drop=FALSE] %*% (sqrt(X$d[p])* t(svdX$u[,p,drop=FALSE]))} 
    return(result)}
   
   
"invMatSqrt" <- function(R,tol = sqrt(.Machine$double.eps)) {
  svdX <- svd(R)
   p <- svdX$d > max(tol * svdX$d[1],0 )
   if(all(p)){ result <- svdX$v %*% (sqrt(1/svdX$d) * t(svdX$u)) } else {      
       result <- svdX$v[,p,drop=FALSE] %*% (sqrt(1/svdX$d[p])* t(svdX$u[,p,drop=FALSE]))
       warning("Sqrt of Penrose inverse used") } 
   result}

