 
 "congeneric.sim" <- 
function(loads = c(0.8, 0.7, 0.6, 0.5),N = NULL,  err=NULL, short=TRUE) { 
 n <- length(loads) 
loading <- matrix(loads, nrow = n) 
error <- diag(1, nrow = n) 
if (!is.null(err)) {diag(error) <- err} else {
 diag(error) <- sqrt(1 - loading^2) }
pattern <- cbind(loading, error) 
colnames(pattern) <- c("theta", paste("e", seq(1:n), sep = "")) 
rownames(pattern) <- c(paste("V", seq(1:n), sep = "")) 
model <- pattern %*% t(pattern)
if(!is.null(N)) {latent <- matrix(rnorm(N * (n + 1)), ncol = (n + 1)) 
       observed <- latent %*% t(pattern) 
       colnames(latent) <-  c("theta", paste("e", seq(1:n), sep = ""))
       if(short) model <- cor(observed) }
if (short) {return(model)}  else {result <- list(model=model,pattern=pattern,r=cor(observed),latent=latent,observed=observed,N=N)
class(result) <- c("psych","sim")
 return(result)} 
 }
 
 
 
 "sim.congeneric" <- 
function(loads = c(0.8, 0.7, 0.6, 0.5),N = NULL,  err=NULL, short=TRUE,categorical=FALSE, low=-3,high=3,cuts=NULL, minor=FALSE,fsmall = c(-.2,.2)) { 
 n <- length(loads) 
 fs <- 0
loading <- matrix(loads, nrow = n) 
if(minor) {
     fsmall  <- c(fsmall,rep(0,n/4))
    fs <- matrix(sample(fsmall,n*floor(n/2),replace=TRUE),ncol=floor(n/2))  
loading <- cbind(loading,fs)}
error <- diag(1, nrow = n ) 
if (!is.null(err)) {diag(error) <- err} else {

 diag(error) <- sqrt(1 - rowSums(loading^2)) }
pattern <- cbind(loading, error) 
if(minor) { colnames(pattern) <- c("theta",paste0("minor",1:NCOL(fs)), paste0("e", seq(1:n)) ) } else {colnames(pattern) <- c("theta", paste0("e", seq(1:n))) }
rownames(pattern) <- c(paste("V", seq(1:n), sep = "")) 
model <- pattern %*% t(pattern)

if(!is.null(N)) {latent <- matrix(rnorm(N * (n + 1 +minor*NCOL(fs))), ncol = (n + 1 +minor*NCOL(fs))) 
       observed <- latent %*% t(pattern) 
       if(is.null(cuts)) {
        if (categorical) {
        
    	observed = round(observed)       #round all items to nearest integer value
		observed[(observed<= low)] <- low     
		observed[(observed>high) ] <- high   
		} } else {
		        temp <- observed
		       	ncuts <- length(cuts)
		 		temp[(observed<= cuts[1])] <- 1 	
		      	if(ncuts > 1) {for (nc in 2:ncuts)  {temp[(observed > cuts[nc-1]) & (observed <= cuts[nc])] <- nc}}
		      	temp[(observed >  cuts[ncuts])] <- ncuts+1 
		      	observed <- temp-1
		}
       colnames(latent) <-  colnames(pattern) #c("theta", paste("e", seq(1:n), sep = ""))
       if(short) model <- cor(observed) }
if (short) {return(model)}  else { if(!is.null(N)) {

result <- list(model=model,pattern=pattern,r=cor(observed),latent=latent,observed=observed,N=N) } else { result<- model} 
class(result) <- c("psych","sim")
 return(result)} 
 }

 