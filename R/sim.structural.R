"sim.structural" <-
function (s=NULL,f=NULL,n=0,mx=NULL,my=NULL,raw=FALSE) {
 require(MASS)
 #first, do the default case of two equally strong predictors
 if(is.null(s)) { 
          s=diag(1,nrow=3)
           s[3,1]  <- s[3,2] <- sqrt(.5)
           s[3,3] <-0}
 
 if(is.null(f) ) {if(!is.null(mx)) { if(is.null(my)) stop("must specify my if mx is specified")
  f <- super.matrix(mx,my)} else {
  f <-matrix(c(               
            .9,0,0,
            .8,0,.0,
            .6,.6,.0,
             0,.8,.0,
             0,-.7,.0,
              0,0,.6,
              0,0,.5,
              0,0,.4),   ncol=3,byrow=TRUE)} }
  
  
 
  if(!is.matrix(f)) f <- as.matrix(f)  #this is the case if doing a congeneric model
  
 
  model <-  f %*% s %*% t(s) %*%  t(f) #the model correlation matrix for oblique factors

  diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(f)[1]
  colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
  if(n>0) {
    mu <- rep(0,nvar)
  	data <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  model <- cor(data) } 
  	reliability <- diag(f %*% t(f))
  if (!raw) {results <- list( model=model,reliability=reliability )} else {
             results <- list( model=model,reliability=reliability,data= data) }
 return(results)}
 
 

 