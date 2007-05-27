"eigen.loadings" <-
function (x) { 
 return(x$vectors %*% sqrt(diag(x$values))) }
 #convert eigen vectors to  principal component loadings by unnormalizing them
 #used if using princomp or princ, not needed for principal
