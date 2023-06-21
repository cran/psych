# A function to create a correlation matrix with a hierarchical structure
# The default values match those of Jensen and Weng
# dropped the call to mvrnorm Nov 28, 2014
# and added back the mu parameter in Sept,2016 (thanks to Alan Robinson for noticing this)
# "sim.hierarchical" <-
# function (gload=NULL,fload=NULL,n=0,raw=FALSE,mu = NULL) {
# cl <- match.call()
#  if(is.null(gload)) gload=matrix(c(.9,.8,.7),nrow=3)
#  if(is.null(fload)) {fload <-matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)}
# 
#   fcor <- gload %*% t(gload)           #the factor correlation matrix
#   diag(fcor) <-1                       #put ones on the diagonal
#   model <-  fload %*% fcor %*% t(fload) #the model correlation matrix for oblique factors
#   diag(model)<- 1                       # put ones along the diagonal
#   nvar <- dim(fload)[1]
#   colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
#   if(n>0) {
#     if(is.null(mu)) mu <- rep(0,nvar)
#   	observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
#   	 the next 3 lines replaces mvrnorm (adapted from mvrnorm, but without the checks)
#                                       eX <- eigen(model)
#                                       observed <-matrix(rnorm(nvar * n),n)
#                                       observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed))
#                                       observed <- t(t(observed) + mu)
#                                       colnames(observed) <- paste("V",1:nvar,sep="")
#   	 r <- cor(observed)
#   	 if(!raw) { 
#   	        result <- list(model=model, r= r, N=n,Call=cl) } else {
#   	        result <- list(model=model, r= r,observed=observed, N=n,Call=cl)}
#   	   class(result) <- c("psych","sim")
#        return(result)
#   	  } else  {return( model) } 
#   
#   }

#An alternative model is to simulate independent factors with cross loadings
#simulate the Thomson Bond'd model to produce a positive manifold without a g factor
#Developed Sept 14, 2019

"sim.bonds" <- function(nvar=9,loads=c(0,0,.5,.6),validity=.8) {
nf <- length(loads)
f <- matrix(0,nrow=(nvar+ nf),ncol=nf)
for (i in 1:nvar) {
   f[i,] <- sample(loads,nf)
   }
    k <- 1
for(i in (nvar+1):NROW(f)) { 
  f[i,k] <- validity
  k <- k+1 }
  colnames(f) <- paste0("F",1:nf)
  rownames(f) <- paste0("V",1:NROW(f)) #fill them all in first
  rownames(f)[(nvar+1) : (nvar+nf)] <- paste0("F",1:nf)
 R <- f %*% t(f)
 diag(R) <- 1
 return(list(R=R,model=f) )
}





# A function to create a correlation matrix with a hierarchical structure
# The default values match those of Jensen and Weng
#completely rewritten, March 15, 2021 to add in the ability to report factor sccores
#modified May 16 to be able to return "items" rather than continuous scales

"sim.hierarchical" <-
function (gload=NULL,fload=NULL,n=0, raw=TRUE, mu = NULL,categorical=FALSE, low=-3,high=3,threshold=NULL) {

cl <- match.call()
#first, prepare the jensen defaults
 if(is.null(gload)) gload=matrix(c(.9,.8,.7),nrow=3) 
 if(is.null(fload)) fload <-matrix(c(.8,.7,.6,rep(0,9),.7,.6,.5,rep(0,9),.6,.5,.4),   ncol=3)
  if(is.null(colnames(fload))) cn <- paste0("F",1:NCOL(fload))
  if(is.null(rownames(fload)))  rn <- paste0("V",1:NROW(fload))
  

  nfact <- NCOL(fload)
  fcor <- gload %*% t(gload)             #the factor correlation matrix
  diag(fcor) <- 1                        #put ones on the diagonal of the factor correlations
  model <-  fload %*% fcor %*% t(fload) #the model correlation matrix for oblique factors
  U <-   diag(sqrt(1-(diag(model))))
  diag(model)<- 1                       # put ones along the diagonal
 
    Ig <- diag(drop(gload))   
    primeload <- fload %*% Ig    #the general factor loadings on each group factor
    fstar  <- sqrt(fload^2 - primeload^2)  #The Schmid Leiman loadings (fload with g removed)
    gstar <- fload %*% gload         #combine them
    gfstar <- cbind(gstar,fstar)   #the factor model with a g factor and orthogonal sl factors
   #house keeping  -- put i some nice labels
    colnames(gfstar) <- c("g",paste0("F",1:nfact,"*"))  #label the lower level as F* to match omega
     nvar <- nrow(fload)
    colnames(model) <- rownames(model)<- paste0("V",1:nvar) 
    n.groups <- NCOL(gfstar )
    if(!is.null(threshold)) {if (length(threshold) < nvar) threshold <- sample(threshold, nvar, replace=TRUE)}
    
  if(n>0) { #make up the data
    
    
    if(is.null(mu)) mu <- rep(0,n.groups)  #don't forget to add this in later
  
     g <-   t(t(matrix(rnorm(n *(n.groups ) ),ncol=n.groups))+mu)  #the factor scores+ factor means
      colnames(g)<- c("g",paste0("F",1:nfact,"*"))
   
    #now generate the lower order variables
     theta <- t(gfstar %*% t(g))  #g are the orthogonal factor scores, theta are the noise free observed
     error <- matrix(rnorm(n*nvar),ncol=nvar)
     observed <-  theta + error %*% (U)  #weight the error by uniqueness of the variables
     
     if(categorical) {#convert from continuous to categorical
         if(is.null(threshold)) {
    	observed = round(observed)       #round all items to nearest integer value
		observed[(observed <= low)] <- low     
		observed[(observed > high) ] <- high 
		} else {
		 i <- 1:nvar
		  observed <- t(t(observed [,i]) > threshold[i]) + 0 }  
		}

     colnames(observed) <- paste0("V",1:nvar)
    
  	 r <- cor(observed)
  	 colnames(r) <- rownames(model)<- paste0("V",1:nvar)
  	 if(!raw) { 
  	        result <- list(model=model, r= r, N=n,Call=cl) } else {
  	        result <- list(model=model, r= r,observed=observed,theta=g, sl=gfstar, N=n,Call=cl)} #include the sl solution
  	   class(result) <- c("psych","sim")
       return(result)
       
  	 } else  {return( model) } #just the modeled matrix
  
  }
  
  