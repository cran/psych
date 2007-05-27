"congeneric.sim" <- 
function(N = 1000, loads = c(0.8, 0.7, 0.6, 0.5), err=NULL, short=TRUE) { 
 n <- length(loads) 
loading <- matrix(loads, nrow = n) 
error <- diag(1, nrow = n) 
if (!is.null(err)) {diag(error) <- err} else {
 diag(error) <- sqrt(1 - loading^2) }
pattern <- cbind(loading, error) 
colnames(pattern) <- c("theta", paste("e", seq(1:n), sep = "")) 
rownames(pattern) <- c(paste("V", seq(1:n), sep = "")) 
latent <- matrix(rnorm(N * (n + 1)), ncol = (n + 1)) 
observed <- latent %*% t(pattern) 
if (short) {return(observed)}  else {result <- list(observed=observed,pattern=pattern)
 return(result)} 
 }