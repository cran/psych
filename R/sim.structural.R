"sim.structure" <- "sim.structural" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,uniq=NULL,raw=TRUE, items = FALSE, low=-2,high=2,d=NULL,cat=5,mu=0) {
 cl <- match.call()
 if(is.null(f)) { if(is.null(fy)) {f <- fx} else {
    f <- superMatrix(fx,fy)} }
  f <- as.matrix(f)
  if(!is.null(Phi)) {if(length(Phi)==1) Phi <- matrix(c(1,Phi,Phi,1),2,2)}
   #these are parameters for simulating items
   nf <- ncol(f)
   nvar <- nrow(f)
if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar/nf-1))
            d <- rep(d,nf)} else {if(length(d)==1) d <- rep(d,nvar)}
	a <- rep(1,nvar)         
  
  if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
                    Phi <- 1}
  if(!is.null(Phi)) {
  model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
} else { model <- f%*% t(f)}
if(is.null(uniq)) {diag(model) <- 1 } else { diag(model) <- uniq  + diag(model)}                      # put ones along the diagonal unless uniq is specified 
  nvar <- dim(f)[1]
  if(is.null(rownames(model))) {colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")} #else {colnames(model) <- rownames(model) <- rownames(fx)}
 
  if(n>0) {
    mu <- rep(mu,nvar)
  	#observed <- mvrnorm(n = n, mu, Sigma=model, tol = 1e-6, empirical = FALSE)
  	 eX <- eigen(model)
                                      observed <- matrix(rnorm(nvar * n),n)
                                      observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed) + mu) 
  	if(items) {observedp <- matrix(t(pnorm(a*t(observed)- d)),n,nvar) 
  	         observed[] <- rbinom(n*nvar, cat, observedp)}
  	  colnames(observed) <- colnames(model)
  r <- cor(observed) 
  } 

  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}
 
  
"sim" <-
function (fx=NULL,Phi=NULL,fy=NULL,alpha=.8,lambda = 0,n=0,mu=NULL,raw=TRUE) {
 cl <- match.call()
##set up some default values 
if(is.null(fx)) {fx <- matrix(c(rep(c(.8,.7,.6,rep(0,12)),3),.8,.7,.6),ncol=4)
   if(is.null(Phi)) {Phi <- diag(1,4,4)
                     Phi <- alpha^abs(row(Phi) -col(Phi)) + lambda^2
                     diag(Phi) <- max((alpha + lambda),1)
	                Phi <- cov2cor(Phi)}
   if(is.null(mu)) {mu <- c(0,.5,1,2)}                   }
    if(is.null(fy)) {f <- fx} else {
    f <- superMatrix(fx,fy)} 
  
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
    
  #	observed <- mvrnorm(n = n, means, Sigma=model, tol = 1e-6, empirical = FALSE)
  	 eX <- eigen(model)
     observed <- matrix(rnorm(nvar * n),n)
     observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed) + rep(means,n))
  r <- cor(observed) } 
  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}
 
 	
"sim.simplex" <-
	function(nvar =12, alpha=.8,lambda=0,beta=1,mu=NULL, n=0) {
	 cl <- match.call()
	R <- matrix(0,nvar,nvar)
	R[] <- alpha^abs(col(R)-row(R))*beta + lambda^2
	diag(R) <- max((alpha * beta) + lambda,1)
	R <- cov2cor(R)
	colnames(R) <- rownames(R) <- paste("V",1:nvar,sep="")
	#require(MASS)
	if(is.null(mu)) {mu <- rep(0,nvar)} 
	if(n>0) {
	#observed.scores <- mvrnorm(n = n, mu, Sigma=R, tol = 1e-6, empirical = FALSE)
	observed <- matrix(rnorm(nvar*n),n)
	 	 eX <- eigen(R)
                observed.scores <- matrix(rnorm(nvar * n),n)
                observed.scores <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed)+mu)
	observed <- cor(observed.scores)
	results <- list(model=R,r=observed,observed=observed.scores)
	results$Call <- cl
	class(results) <- c("psych", "sim")} else {results <- R}
	  
	results
	}
	

#simulate major and minor factors
"sim.minor" <-
function(nvar=12,nfact=3,n=0,g=NULL,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE) {
if(is.null(fbig)) {loads <- c(.8,.6) } else {loads <- fbig}

loads <- sample(loads,nvar/nfact,replace=TRUE)
if(nfact == 1) {fx <- matrix(loads,ncol=1)} else {fx <- matrix(c(rep(c(loads,rep(0,nvar)),(nfact-1)),loads),ncol=nfact)}
if(bipolar) fx <- 2*((sample(2,nvar,replace=TRUE) %%2)-.5) * fx
if(!is.null(g)) {if (length(g) < nvar) {g <- sample(g,nvar,replace=TRUE)}
        fx <- cbind(g,fx)
        }
fsmall  <- c(fsmall,rep(0,nvar/4))
fs <- matrix(sample(fsmall,nvar*floor(nvar/2),replace=TRUE),ncol=floor(nvar/2))  
fload <- cbind(fx,fs)
if(is.null(g)) {
colnames(fload) <- c(paste("F",1:nfact,sep=""),paste("m",1:(nvar/2),sep=""))} else {
    colnames(fload) <- c("g",paste("F",1:nfact,sep=""),paste("m",1:(nvar/2),sep=""))}
rownames(fload) <- paste("V",1:nvar,sep="")
results <- sim(fload,n=n)
         results$fload <- fload
 class(results) <- c("psych", "sim")
 return(results)
}


#simulate various structures and summarize them
"sim.omega" <-
function(nvar=12,nfact=3,n=500,g=NULL,sem=FALSE,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE,om.fact=3,flip=TRUE,option="equal",ntrials=10) {
results <- matrix(NaN,nrow=ntrials,ncol=12)
colnames(results) <- c("n","om.model","omega","ev.N","e.f1","omega.f1","Beta","omegaCFA","omegaSem","rms","RMSEA","coeff.v")
for (i in 1:ntrials) {
x <- try(sim.minor(nvar=nvar,nfact=nfact,n=n,g=g,fbig=fbig,fsmall=fsmall,bipolar=bipolar))
if(is.null(g)) {omega.model <- 0} else {gsum <- colSums(x$fload)[1]
    omega.model <- gsum^2/sum(x$model)}
results[i,"om.model"] <- omega.model
observed.cor <- cor(x$observed)
ev <- eigen(observed.cor)$values
f1 <- fa(observed.cor)$loadings
om.fa <- sum(f1)^2/sum(observed.cor)
e.f1 <- sum(f1^2)/nvar
    sem.model <- omega.sem(x$fload,sl=TRUE,nf=nfact)  #this is the model based upon the true values
    if(sem) {if(!requireNamespace('sem')) {stop("You must have the sem package installed to use omegaSem")} else {sem.om <- try(sem(model=sem.model,S=observed.cor, N=n))} 
 
  omega.cfa <- omegaFromSem(observed.cor,sem.om,flip=flip)
    if(omega.cfa$omega >1) omega.cfa$omega <- NA
  results[i,"omegaCFA"] <- omega.cfa$omega
  } else {omega.cfa <- NULL}
results[i,"n"] <- n
if(n > 0) {
   if (sem) {om <- try(omegaSem(x$observed,om.fact,flip=flip,plot=FALSE,option=option))} else {
om <- try(omega(x$observed,om.fact,flip=flip,plot=FALSE,option=option))}
          ic <- suppressWarnings(ICLUST(x$observed,1,plot=FALSE))} else {
           if (sem) {om <- try(omegaSem(x$model,om.fact,flip=flip,plot=FALSE,option=option))} else {om <- try(omega(x$model,om.fact,flip=flip,plot=FALSE,option=option))
           if(class(om)=="try-error") {message("Error in sem. iteration = ",i)
            om <- NA
            next}}
          ic <- suppressWarnings(ICLUST(x$model,1,plot=FALSE))}
results
if(sem) {results[i,"omega"] <- om$omegaSem$omega_h
                loads <- om$omegaSem$schmid$sl
     
} else {results[i,"omega"] <- om$omega_h
loads <- om$schmid$sl
}  
 p2 <- loads[,ncol(loads)]
         	mp2 <- mean(p2)
         	vp2 <- var(p2)       	
       
#results[i,"p2"] <-  mp2
#results[i,"p2.sd"] <-  sqrt(vp2)

results[i,"coeff.v"] <- sqrt(vp2)/mp2
results[i,"Beta"] <-  ic$beta
results[i,"ev.N"] <- ev[1]/nvar
results[i,"e.f1"] <- e.f1
results[i,"omega.f1"] <- om.fa


if(sem) {
        if(!is.null(om$omegaSem$schmid$RMSEA)) {results[i,"RMSEA"] <-om$omegaSem$schmid$RMSEA[1]} else {results[i,"RMSEA"] <- NA}
         if(!is.null(om$omegaSem$schmid$rms))results[i,"rms"] <- om$omegaSem$schmid$rms
        results[i,"omegaSem"] <- om$omega.efa$omega
        if(results[i,"omegaSem"] > 1) {warning("Bad result from sem   case = ",i) 
                             results[i,"omegaSem"] <- NA}
        } else { 
if(!is.null(om$schmid$RMSEA)) {results[i,"RMSEA"] <- om$schmid$RMSEA[1]} else {results[i,"RMSEA"] <- NA}
if(!is.null(om$schmid$rms)) results[i,"rms"] <- om$schmid$rms
     results[i,"omegaSem"] <- NA}

}

if(n==0) {results <- results[,-which(colnames(results)=="RMSEA")] #drop RMSEA if there are no cases
if(!sem) results <- results[,-which(colnames(results)=="omegaSem")] 
} else {if(!sem) results <- results[,-which(colnames(results)=="omegaSem")] }

return(results)
}

"sim.omega.2" <- function(nvar=12,nfact=3,n=c(100,200,400,800),g=c(0,.1,.2,.3,.4,.5),sem=TRUE,fbig=c(.7,.6),fsmall=c(-.2,.2),bipolar=FALSE,om.fact=3,ntrials=10) {
result <- list() 
k <- 1
progressBar(k,length(n)*length(g),"sim.omega.2")
for (ni in 1:length(n)) {
for (gi in 1:length(g)) {
result[[k]] <- sim.omega(nvar=nvar,nfact=nfact,n=n[ni],g =g[gi],fbig=fbig,fsmall=fsmall,bipolar=bipolar,ntrials=ntrials,om.fact=om.fact,sem=sem)
k <- k+1
}
}
cnames <- colnames(result[[1]])
#result <- unlist(result)
#if(sem) {result <- matrix(result,ncol=10,byrow=TRUE)} else {result <- matrix(result,ncol=9,byrow=TRUE) }
#colnames(result) <- cnames

return(result)
}


"sim.general" <- 
function(nvar=9,nfact=3, g=.3,r=.3,n=0) {
#require(MASS)
  r1 <- matrix(r,nvar/nfact,nvar/nfact)
  R <- matrix(g,nvar,nvar)
  rf <- superMatrix(r1,r1)
  if(nfact>2) {for (f in 1:(nfact-2)){
  rf <- superMatrix(r1,rf)}}
  R <- R + rf
  diag(R) <- 1
  colnames(R) <- rownames(R) <- paste((paste("V",1:(nvar/nfact),sep="")),rep(1:nfact,each=(nvar/nfact)),sep="gr")
  if(n > 0) {#x <-  mvrnorm(n = n, mu=rep(0,nvar), Sigma = R, tol = 1e-06,empirical = FALSE)
   eX <- eigen(R)
                                      x <- matrix(rnorm(nvar * n),n)
                                      x <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(x))
   return(x)} else {
  return(R)} }