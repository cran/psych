"sim.congeneric" <- 
function(N = 1000, loads = c(0.8, 0.7, 0.6, 0.5), err=NULL, short=TRUE) { 
 n <- length(loads) 
loading <- matrix(loads, nrow = n) 
error <- diag(1, nrow = n) 
if (!is.null(err)) {diag(error) <- err} else {
 diag(error) <- sqrt(1 - loading^2) }
pattern <- cbind(loading, error) 
colnames(pattern) <- c("theta", paste("e", seq(1:n), sep = "")) 
rownames(pattern) <- c(paste("V", seq(1:n), sep = "")) 
model <- pattern %*% t(pattern)
latent <- matrix(rnorm(N * (n + 1)), ncol = (n + 1)) 
observed <- latent %*% t(pattern) 
colnames(latent) <-  c("theta", paste("e", seq(1:n), sep = "")) 
if (short) {return(model)}  else {result <- list(model=model,pattern=pattern,observed=observed,latent=latent)
 return(result)} 
 }