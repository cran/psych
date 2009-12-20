"make.structural" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,raw=FALSE) {
 require(MASS)
 #first, do the default case of two equally strong predictors
 if(is.null(Phi)) { 
          Phi=diag(1,nrow=3)
           Phi[3,1]  <- Phi[3,2] <- sqrt(.5)
           Phi[3,3] <-0}
 
 if(is.null(f) ) {if(!is.null(fx)) { if(is.null(fy)) stop("must specify my if mx is specified")
  f <- super.matrix(fx,fy) } else {   #create a default example
  fx <-matrix(c( .9,.8,.6,rep(0,4),.6,.8,-.7),ncol=2)              
  fy <- c(.6,.5,.4)  
  f <- super.matrix(fx,fy)
            }
            }
  #this works for asymetric Phi matrices, which are regression weights
  if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
                    Phi <- 1}
  
if(!isSymmetric(Phi)) {  model <-  f %*% Phi %*% t(Phi) %*%  t(f)} else { 
                          model <-   f %*% Phi  %*%  t(f)} #the model correlation matrix for oblique factors

  diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(f)[1]
  if(is.null(colnames(model))) {colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")}
  if(n>0) {
    mu <- rep(0,nvar)
  	data <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  r <- cor(data) } 
  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= data,N=n) } }
  class(results) <- c("psych", "sim")
 return(results)}
 
 

 