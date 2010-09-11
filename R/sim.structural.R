"sim.structure" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,raw=FALSE) {
 cl <- match.call()
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
  if(is.null(rownames(model))) {colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")} #else {colnames(model) <- rownames(model) <- rownames(fx)}
 
  if(n>0) {
    mu <- rep(0,nvar)
  	observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  r <- cor(observed) } 
  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}
 
 
"sim.structural" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,raw=FALSE) {
require(MASS)
 cl <- match.call()
 if(is.null(f)) { if(is.null(fy)) {f <- fx} else {
    f <- super.matrix(fx,fy)} }
  
  if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
                    Phi <- 1}
  if(!is.null(Phi)) {
  
  model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
} else { model <- f%*% t(f)}
diag(model)<- 1                       # put ones along the diagonal
  nvar <- dim(f)[1]
  if(is.null(rownames(model))) {colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")} #else {colnames(model) <- rownames(model) <- rownames(fx)}
  if(n>0) {
    mu <- rep(0,nvar)
  	observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  r <- cor(observed) } 
  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}
 
 
"sim" <-
function (fx=NULL,Phi=NULL,fy=NULL,n=0,mu=NULL,raw=TRUE) {
 cl <- match.call()
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
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}
 
 	
	"sim.simplex" <-
	function(nvar =12, r=.8,mu=NULL, n=0) {
	 cl <- match.call()
	R <- matrix(0,nvar,nvar)
	R[] <- r^abs(col(R)-row(R))
	colnames(R) <- rownames(R) <- paste("V",1:nvar,sep="")
	require(MASS)
	if(is.null(mu)) {mu <- rep(0,nvar)} 
	if(n>0) {
	observed.scores <- mvrnorm(n = n, mu, Sigma=R, tol = 1e-6, empirical = FALSE)
	observed <- cor(observed.scores)
	results <- list(model=R,r=observed,observed=observed.scores)
	results$Call <- cl
	class(results) <- c("psych", "sim")} else {results <- R}
	  
	results
	}
	

#simulate major and minor factors
"sim.minor" <-
function(nvar=12,nfact=3,n=0,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE) {
if(is.null(fbig)) {loads <- c(.8,.6) } else {loads <- fbig}
loads <- sample(loads,nvar/nfact,replace=TRUE)
if(nfact == 1) {fx <- matrix(loads,ncol=1)} else {fx <- matrix(c(rep(c(loads,rep(0,nvar)),(nfact-1)),loads),ncol=nfact)}
if(bipolar) fx <- 2*((sample(2,nvar,replace=TRUE) %%2)-.5) * fx
fsmall  <- c(fsmall,rep(0,nvar/4))
fs <- matrix(sample(fsmall,nvar*floor(nvar/2),replace=TRUE),ncol=floor(nvar/2))  
fload <- cbind(fx,fs)
colnames(fload) <- c(paste("F",1:nfact,sep=""),paste("m",1:(nvar/2),sep=""))
rownames(fload) <- paste("V",1:nvar,sep="")
results <- sim(fload,n=n)
         results$fload <- fload
 class(results) <- c("psych", "sim")
 return(results)
}


#similuate various structures and summarize them
"sim.omega" <-
function(nvar=12,nfact=3,n=0,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE,om.fact=3,flip=TRUE,option="equal",ntrials=10) {
results <- matrix(NA,nrow=ntrials,ncol=7)
for (i in 1:ntrials) {
x <- sim.minor(nvar=nvar,nfact=nfact,n=n,fbig=fbig,fsmall=fsmall,bipolar=bipolar)
if(n > 0) {
om <- try(omega(x$observed,om.fact,flip=flip,plot=FALSE,option=option))
          ic <- ICLUST(x$observed,1,plot=FALSE)} else {om <- try(omega(x$model,om.fact,flip=flip,plot=FALSE,option=option))
          ic <- ICLUST(x$model,1,plot=FALSE)}

results[i,1] <- om$omega_h
     loads <- om$schmid$sl
   p2 <- loads[,ncol(loads)]
         	mp2 <- mean(p2)
         	vp2 <- var(p2)
         	
       
results[i,2] <-  mp2
results[i,3] <-  sqrt(vp2)
results[i,4] <- sqrt(vp2)/mp2
results[i,5] <-ic$beta
if(!is.null(om$schmid$RMSEA)) {results[i,6] <-om$schmid$RMSEA[1]} else {results[i,6] <- NA}
if(!is.null(om$schmid$rms))results[i,7] <- om$schmid$rms


}
colnames(results) <- c("omega","mean p2","sd p2","coeff v","Beta","RMSEA","rms")
return(results)
}
