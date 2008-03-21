"eigen.loadings" <-
function (x) { 
if(!is.null(x$vector)) {
ans <- x$vectors %*% sqrt(diag(x$values))
colnames(ans) <- rownames(ans) <- rownames(x$vector)
 return(ans)
 } else 
 if(!is.null(x$loadings)) {
 ans <- x$loadings %*% diag(x$sdev)
 rownames(ans) <- rownames(x$loadings)
 colnames(ans) <- colnames(x$loadings)
 return(ans)
 }
 }
 #convert eigen vectors to  principal component loadings by unnormalizing them
 
