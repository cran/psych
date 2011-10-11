"kaiser" <- function(f,rotate="oblimin") {
 if((!is.matrix(f)) && (!is.data.frame(f)))  {f <- as.matrix(f$loadings)} else {f <- as.matrix(f)}

h2 <- diag(f %*% t(f))
weighted <- f/sqrt(h2)
rotated <- do.call(rotate,list(weighted))
normalized <- rotated$loadings * sqrt(h2)
rotated$loadings <- normalized
return(rotated)}
