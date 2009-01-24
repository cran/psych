# A function to create a correlation matrix with a hierarchical structure
"sim.hierarchical" <-
function (gload=NULL,fload=NULL,n=0,raw=FALSE) {
 require(MASS)
 if(is.null(gload)) gload=matrix(c(.9,.8,.7),nrow=3)
 if(is.null(fload)) {fload <-matrix(c(.9,.8,.7,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)}

  fcor <- gload %*% t(gload)           #the factor correlation matrix
  diag(fcor) <-1                       #put ones on the diagonal
  model <-  fload %*% fcor %*% t(fload) #the model correlation matrix for oblique factors
  diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(fload)[1]
  colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
  if(n>0) {
    mu <- rep(0,nvar)
  	observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  	 r <- cor(observed)
  	 if(!raw) { 
  	        result <- list(model=model, r= r, N=n) } else {
  	        result <-list( model=model, r= r,observed=observed, N=n)}
  	   class(result) <- c("psych","sim")
       return(result)
  	  } else  {return( model) } 
  
  }
