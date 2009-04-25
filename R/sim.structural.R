"sim.structure" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,raw=FALSE) {
require(MASS)

 if(is.null(f)) { if(is.null(fy)) {f <- fx} else {
    f <- super.matrix(fx,fy)} }
  
            
  
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
 
 
"sim.structural" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,raw=FALSE) {
require(MASS)

 if(is.null(f)) { if(is.null(fy)) {f <- fx} else {
    f <- super.matrix(fx,fy)} }
  
            
  
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
 
 
"sim" <-
function (fx=NULL,Phi=NULL,fy=NULL,n=0,mu=NULL,raw=FALSE) {
require(MASS)
#set up some default values 
if(is.null(fx)) {fx <- matrix(c(rep(c(.8,.7,.6,rep(0,12)),3),.8,.7,.6),ncol=4)
   if(is.null(Phi)) {Phi <- diag(1,4,4)
                     Phi <- .8^abs(row(Phi) -col(Phi))}
   if(is.null(mu)) {mu <- c(0,.5,1,2)}                   }
    if(is.null(fy)) {f <- fx} else {
    f <- super.matrix(fx,fy)} 
  
   if(is.null(mu)) {mu <- rep(0,ncol(fx))} 
   
     means <- fx %*% mu        
  
  if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
                    Phi <- 1}
  if(!is.null(Phi)) {
  model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
} else { model <- f%*% t(f)}
diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(f)[1]
  colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
  if(n>0) {
    
  	observed <- mvrnorm(n = n, means, Sigma=model, tol = 1e-6, empirical = FALSE)
  r <- cor(observed) } 
  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
  class(results) <- c("psych", "sim")
 return(results)}
 
 	
	"sim.simplex" <-
	function(nvar =12, r=.8,mu=NULL, n=0) {
	R <- matrix(0,nvar,nvar)
	R[] <- r^abs(col(R)-row(R))
	require(MASS)
	if(is.null(mu)) {mu <- rep(0,nvar)} 
	if(n>0) {
	observed.scores <- mvrnorm(n = n, mu, Sigma=R, tol = 1e-6, empirical = FALSE)
	observed <- cor(observed.scores)
	results <- list(model=R,r=observed,observed=observed.scores) } else {results <- R}
	results
	}