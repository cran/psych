# A function to create a correlation matrix with a hierarchical structure
# The default values match those of Jensen and Weng
#dropped the call to mvrnorm Nov 28, 2014
"sim.hierarchical" <-
function (gload=NULL,fload=NULL,n=0,raw=FALSE,mu = NULL) {
cl <- match.call()
# require(MASS)
 if(is.null(gload)) gload=matrix(c(.9,.8,.7),nrow=3)
 if(is.null(fload)) {fload <-matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)}

  fcor <- gload %*% t(gload)           #the factor correlation matrix
  diag(fcor) <-1                       #put ones on the diagonal
  model <-  fload %*% fcor %*% t(fload) #the model correlation matrix for oblique factors
  diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(fload)[1]
  colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
  if(n>0) {
   # if(is.null(mu)) mu <- rep(0,nvar)
  	#observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  	 #the next 3 lines replaces mvrnorm (adapted from mvrnorm, but without the checks)
                                      eX <- eigen(model)
                                      observed <- matrix(rnorm(nvar * n),n)
                                      observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed))
                                      colnames(observed) <- paste("V",1:nvar,sep="")
  	 r <- cor(observed)
  	 if(!raw) { 
  	        result <- list(model=model, r= r, N=n,Call=cl) } else {
  	        result <- list(model=model, r= r,observed=observed, N=n,Call=cl)}
  	   class(result) <- c("psych","sim")
       return(result)
  	  } else  {return( model) } 
  
  }


