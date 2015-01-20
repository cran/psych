"kaiser" <- function(f,rotate="oblimin") {
 if((!is.matrix(f)) && (!is.data.frame(f)))  {f <- as.matrix(f$loadings)} else {f <- as.matrix(f)}
if(!requireNamespace('GPArotation')) stop('GPArotation is required for the Kaiser normalization')
h2 <- diag(f %*% t(f))
weighted <- f/sqrt(h2)
#rotated <- call(paste('GPArotation',rotate,sep="::"),list(weighted))
rotated <- do.call(getFromNamespace(rotate,'GPArotation'),list(weighted))
normalized <- rotated$loadings * sqrt(h2)
rotated$loadings <- normalized
class(rotated) <- c("psych","fa")
return(rotated)}
