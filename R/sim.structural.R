"sim.structural" <-
function (fx=NULL,fy=NULL,Phi=NULL,f=NULL,n=0,raw=FALSE) {
require(MASS)

 if(is.null(f)) {if(!is.null(fx)) { if(is.null(fy)) stop("must specify ymodel if xmodel is specified")
    f <- super.matrix(fx,fy) 
    }}
            
  
  if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
                    Phi <- 1}
  if(!is.null(Phi)) {
  model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
} else { model <- f%*% t(f)}
diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(f)[1]
  colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
  if(n>0) {
    mu <- rep(0,nvar)
  	observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  r <- cor(observed) } 
  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
  class(results) <- c("psych", "sim")
 return(results)}
 
 

 