"sim.structure" <- "sim.structural" <-
function (fx=NULL,Phi=NULL,fy=NULL,f=NULL,n=0,uniq=NULL,raw=TRUE, items = FALSE, low=-2,high=2,d=NULL,cat=5,mu=0) {
 cl <- match.call()
 if(is.null(fx)) {fx <- matrix(c(rep(c(.8,.7,.6,rep(0,12)),3),.8,.7,.6),ncol=4)
             colnames(fx) <- paste0("F",1:4)
             rownames(fx) <- paste0("V",1:12)
            }  #added default 10/21/23
            
 if(is.null(f)) { if(is.null(fy)) {f <- fx} else {
    f <- superMatrix(fx,fy)} }
  f <- as.matrix(f)
  if(!is.null(Phi)) {if(length(Phi)==1) Phi <- matrix(c(1,Phi,Phi,1),2,2)}
   #these are parameters for simulating items
   nf <- ncol(f)
   nvar <- nrow(f)
#if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar/(nf-1)))
         #   d <- rep(d,nvar)} else {if(length(d)==1) d <- rep(d,nvar)}
	#a <- rep(1,nvar)         
  
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
                 theta <- matrix(rnorm(nvar * n),n)   #random values
                observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(theta) + mu)   #this way theta and observed not identical
                 theta <- observed   #but seems to make them identical
  	if(items) { #the next 2 lines were dropped in 2021 because it was adding noise twice
  	         #observedp <- matrix(t(pnorm(a*t(observed)- d)),n,nvar)    #this is not useful
  	         #observed[] <- rbinom(n*nvar, cat, observedp) #this puts in error again
  	         #put in the low and high values again --dropped in 2021
  	         observed[observed <low] <- low
  	         observed[observed > high] <- high
  	         range.ob <- range(observed)
  	         observed <- (observed - range.ob[1])/(range.ob[2]- range.ob[1]) 
  	         observed<- round((cat-1) * observed) + 1}
  	  colnames(observed) <- colnames(model)
  r <- cor(observed) 
  } 

  	reliability <- diag(f %*% t(f))
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,theta=theta, N=n) } }
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}
 
 # Replaced 1/8/21 with a version that reports theta 
# "sim" <-
# function (fx=NULL,Phi=NULL,fy=NULL,alpha=.8,lambda = 0,n=0,mu=NULL,raw=TRUE) {
#  cl <- match.call()
# ##set up some default values 
# if(is.null(fx)) {fx <- matrix(c(rep(c(.8,.7,.6,rep(0,12)),3),.8,.7,.6),ncol=4)
#    if(is.null(Phi)) {Phi <- diag(1,4,4)
#                      Phi <- alpha^abs(row(Phi) -col(Phi)) + lambda^2
#                      diag(Phi) <- max((alpha + lambda),1)
# 	                Phi <- cov2cor(Phi)}
#    if(is.null(mu)) {mu <- c(0,.5,1,2)}                   }
#     if(is.null(fy)) {f <- fx} else {
#     f <- superMatrix(fx,fy)} 
#   
#    if(is.null(mu)) {mu <- rep(0,ncol(fx))} 
#    
#      means <- fx %*% mu        
#   
#   if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
#                     Phi <- 1}
#   if(!is.null(Phi)) {
#   model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
# } else { model <- f%*% t(f)}
# diag(model)<- 1                       # put ones along the diagonal
#   nvar <- dim(f)[1]
#   colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")
#   if(n>0) {
#     
#   #	observed <- mvrnorm(n = n, means, Sigma=model, tol = 1e-6, empirical = FALSE)
#   	 eX <- eigen(model)
#      observed <- matrix(rnorm(nvar * n),n)
#      observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed) + rep(means,n))
#   r <- cor(observed) } 
#   	reliability <- diag(f %*% t(f))
#   if(n<1) {results <- list(model=model,reliability=reliability) } else {
#   if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
#              results <- list( model=model,reliability=reliability,r=r,observed= observed,N=n) } }
#   results$Call <- cl
#   class(results) <- c("psych", "sim")
#  return(results)}
#  


#####
#Added 1/7/21  to report the generating theta
#corrected 1/29/21 to correctly report generating theta (with correlations)
####
 	
# "sim" <-
# function (fx=NULL,Phi=NULL,fy=NULL,alpha=.8,lambda = 0,n=0,mu=NULL,raw=TRUE,theta=TRUE) {
#  cl <- match.call()
#  
# ##set up some default values 
# if(is.null(fx)) {fx <- matrix(c(rep(c(.8,.7,.6,rep(0,12)),3),.8,.7,.6),ncol=4)}
#     if(is.vector(fx)) {fx <- as.matrix(fx) } #this is the case if doing a congeneric model
#                    
#   if(is.vector(fy)) {fy <- as.matrix(fy)}
#      
#     nvar <- NROW(fx)
#     nfact <- NCOL(fx)
#    if(is.null(Phi)& is.null(fy)) {Phi <- diag(nfact)
#                      Phi <- alpha^abs(row(Phi) -col(Phi)) + lambda^2
#                      diag(Phi) <- max((alpha + lambda),1)
# 	                Phi <- cov2cor(Phi)}    #put in a simplex structure as default
#                  
#     if(is.null(fy)) {f <- fx} else {
#     f <- superMatrix(fx,fy)}
#     
#      if(is.null(mu))  mu <- as.matrix(0:(NCOL(f)-1),drop=FALSE)
#      if(is.null(Phi)) Phi <- diag(NCOL(f)) 
#   
#    if(is.null(mu)) {mu <- as.matrix(rep(0,NCOL(f)),drop=FALSE)} 
#    
#      means <- f %*% mu        
#   
# 
#      model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
#    
#    
#      nvar <- NROW(f) #we do this here because we have made f a supermatrix of fx and fy   
#      colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")  
#      
#         
#  if((n > 0)) { 
#      nfactors <- NCOL(f)
#      theta <- matrix(rnorm(n * nfactors),n)   #make up some random data 
#      
#      eX <- eigen(model)  #the model is actually the covariance matrix
#       t.scores <-  t(eX$vectors[,1:nfactors] %*% diag(sqrt(pmax(eX$values[1:nfactors], 0)))    %*%  t(theta)  + rep(means,n)) #the true scores for each variable
#       U <- diag(sqrt(1-diag(model)))
#       error <- t(U %*% matrix(rnorm(n * nvar),nrow=nvar))
#       observed <- t.scores + error
#       colnames(observed) <- paste0("V",1:ncol(observed))
#        r <- cor(observed)
#       if(!is.null(Phi)) { exPhi <- eigen(Phi)
#       # f.scores <- f.scores %*% Phi
# 
#         #f.scores <-t((exPhi$vectors %*% diag(sqrt(exPhi$values)) %* t(f.scores)) + mu)}
#         theta <-    t(exPhi$vectors %*% diag(sqrt(exPhi$values)) %*% t(theta) + rep(mu,n) ) 
#      }
#      } #  else {
#      
# 
#   diag(model)<- 1                       # put ones along the diagonal
#   	reliability <- diag(f %*% t(f))
#   if(n < 1) {results <- list(model=model,reliability=reliability) } else {
#   if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
#              results <- list( model=model,reliability=reliability,r=r, observed= observed,theta=theta,Phi=Phi, N=n) } }
#   results$Call <- cl
#   class(results) <- c("psych", "sim")
#  return(results)} 
#  
#  
 
 
 #complete rewrite to allow the true factor scores (theta) to be reported
#rewritten March, 2021
#added threshold option April 2023
#added the items and low, cat options to make a more general simulation function 10/20/23 
 "sim" <-
function (fx=NULL,Phi=NULL,fy=NULL,alpha=.8,lambda = 0,n=0,mu=NULL,raw=TRUE,threshold=NULL,items = FALSE, low=-2,high=2,cat=5) {
 cl <- match.call()
 
##set up some default values 
# 4 correlated factors growing over time
# The four factors have a simplex structure

  
    
if(is.null(fx)) {fx <- matrix(c(rep(c(.8,.7,.6,rep(0,12)),3),.8,.7,.6),ncol=4)
   	if(is.null(Phi)) {Phi <- diag(NCOL(fx))         #create a simplex of the factors
                     Phi <- alpha^abs(row(Phi) -col(Phi)) + lambda^2
                     diag(Phi) <- max((alpha + lambda),1)
	                Phi <- cov2cor(Phi)}
   if(is.null(mu)) {mu <- seq(0,NCOL(fx)-1)} else {if(length(mu)< NCOL(fx)) mu <- sample(mu,NCOL(fx),replace=TRUE)}                  
      }    #end of default
    if(is.null(fy)) {f <- fx} else {
    f <- superMatrix(fx,fy)} 
     nvar <- NROW(f)
     nfactors <- NCOL(f)

    if(is.vector(f)) {f <- as.matrix(f)  #this is the case if doing a congeneric model
                    Phi <- 1}
   if(is.null(mu)) {mu <- rep(0,NCOL(f))} 
   
     means <- f  %*% mu  
        
  
  if(is.null(Phi)) {Phi <-diag(NCOL(f))}
        model <-  f %*% Phi %*%  t(f) #the model correlation matrix for oblique factors
          U <- diag(sqrt(1-diag(model)))   #the uniquensses 
          diag(model)<- 1                  # put ones along the diagonal
   
     colnames(model) <- rownames(model) <- paste("V",1:nvar,sep="")  
     
 #Create factor scores and raw data       
 if((n>0) ) { 
        if(!is.null(threshold)) {if (length(threshold) < nvar) threshold <- sample(threshold, nvar, replace=TRUE)}
     theta <- matrix(rnorm(n * nfactors),ncol=nfactors)   #the orthogonal factors
     et <- eigen(Phi)
  
     theta <- t(et$vectors %*% diag(sqrt(et$values ))%*% t(theta))  #the correlated factors
     observed <- theta %*% t(f)
     colnames(theta) <- colnames(f) 
      error <- t(U %*% matrix(rnorm(n * nvar),nrow=nvar))
      observed <- t(t(observed + error) + rep(means,n))  #observed are factors * loadings + error + means
      	if(items) { 
  	         observed[observed <low] <- low
  	         observed[observed > high] <- high
  	         range.ob <- range(observed)
  	         observed <- (observed - range.ob[1])/(range.ob[2]- range.ob[1]) 
  	         observed<- round((cat-1) * observed) + 1}
      
       if(!is.null(threshold)) {
        i <- 1:nvar
		  observed <- t(t(observed [,i]) > threshold[i]) + 0 } 
       r <- cor(observed)  #will approximate the model
     }  
  	reliability <- diag(f %*% t(f))   
  if(n<1) {results <- list(model=model,reliability=reliability) } else {
  if (!raw) {results <- list( model=model,reliability=reliability,r=r,N=n )} else {
             results <- list( model=model,reliability=reliability,r=r,observed= observed,theta=theta,N=n) } }
  results$Call <- cl
  class(results) <- c("psych", "sim")
 return(results)}

 
 ###########################	
"sim.simplex" <-
	function(nvar =12, alpha=.8,lambda=0,beta=1,mu=NULL, n=0,threshold=NULL) {
	 cl <- match.call()
	R <- matrix(0,nvar,nvar)
	R[] <- alpha^abs(col(R)-row(R))*beta + lambda^2
	diag(R) <- max((alpha * beta) + lambda,1)
	R <- cov2cor(R)
	colnames(R) <- rownames(R) <- paste("V",1:nvar,sep="")
	
	if(is.null(mu)) {mu <- rep(0,nvar)} 
	
	if(n>0) {
	#observed.scores <- mvrnorm(n = n, mu, Sigma=R, tol = 1e-6, empirical = FALSE)
	observed <- matrix(rnorm(nvar *n),n)
	 	 eX <- eigen(R)
                observed.scores <- matrix(rnorm(nvar * n),n)
                observed.scores <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed)+mu)
    if(!is.null(threshold)) {if (length(threshold) < nvar) threshold <- sample(threshold, nvar, replace=TRUE)
             i <- 1:nvar
		  observed.scores <- t(t(observed.scores [,i]) > threshold[i]) + 0 } 
    colnames(observed.scores)<- paste0("V",1:nvar)           
	observed <- cor(observed.scores)
	results <- list(model=R,r=observed,observed=observed.scores) 
	results$Call <- cl
	class(results) <- c("psych", "sim")} else {results <- R}
	  
	results
	}
	

#simulate major and minor factors
#modified October 18, 2022 to allow specification of g and fbig for all variables
#further modified December 4, 2023 to allow specificaton of Ph
"sim.minor" <-
function(nvar=12, nfact=3,n=0,g=NULL,fbig=NULL,fsmall = c(-.2,.2),n.small=NULL,Phi=NULL, bipolar=TRUE, threshold=NULL,
  items = FALSE, low=-2,high=2,cat=5  ) {
if(is.null(fbig)) {loads <- c(.8,.6) 
  loads <- sample(loads,nvar/nfact,replace=TRUE)


if(nfact == 1) {fx <- matrix(loads,ncol=1)} else {fx <- matrix(c(rep(c(loads,rep(0,nvar)),(nfact-1)),loads),ncol=nfact)}
if(bipolar) fx <- 2*((sample(2,nvar,replace=TRUE) %%2)-.5) * fx
} else {fx <- fbig
     nvar <- NROW(fx)
     nfact <- NCOL(fx)}
if(!is.null(g)) {if (length(g) < nvar) {g <- sample(g,nvar,replace=TRUE)}
        fx <- cbind(g,fx)
        }
 
if( !is.null(fsmall)) {
if(is.null(n.small)) n.small=nvar/4
fsmall  <- c(fsmall,rep(0,n.small))
fs <- matrix(sample(fsmall,nvar*floor(n.small),replace=TRUE),ncol=floor((n.small)))  
fload <- cbind(fx,fs)

if(is.null(g)) {
colnames(fload) <- c(paste("F",1:nfact,sep=""),paste("m",1:(n.small),sep=""))} else {
    colnames(fload) <- c("g",paste("F",1:nfact,sep=""),paste("m",1:(n.small),sep=""))}
    } else {
     fload <- fx
      colnames(fload) <- paste0("F",1:nfact)
      }
rownames(fload) <- paste0("V",1:nvar)
sanity.check <- h2  <- apply(fload,1,function(x) sum(x^2))
if(max(sanity.check) >  1.0) {print(cbind(fload, h2))
                       stop("Model is impossible.  Communalities (h2) exceed 1. ")}
if(!is.null(Phi)) {Phi <- superMatrix(Phi,diag(n.small))}  #added 12/3/23
results <- sim(fload,n=n, Phi=Phi,threshold=threshold,  items = items, low=low,high=high,cat=cat )
         results$fload <- fload
 class(results) <- c("psych", "sim")
 return(results)
}


#simulate various structures and summarize them
"sim.omega" <-
function(nvar=12,nfact=3,n=500,g=NULL,sem=FALSE,fbig=NULL,fsmall = c(-.2,.2),bipolar=TRUE,om.fact=3,flip=TRUE,option="equal",ntrials=10,threshold=NULL) {
results <- matrix(NaN,nrow=ntrials,ncol=12)
colnames(results) <- c("n","om.model","omega","ev.N","e.f1","omega.f1","Beta","omegaCFA","omegaSem","rms","RMSEA","coeff.v")

#iterate the entire process n.trials times
for (i in 1:ntrials) {
#make up some data using the Tucker notion of big and small factors
x <- try(sim.minor(nvar=nvar,nfact=nfact,n=n,g=g,fbig=fbig,fsmall=fsmall,bipolar=bipolar,threshold=threshold))
if(is.null(g)) {omega.model <- 0} else {gsum <- colSums(x$fload)[1]
    omega.model <- gsum^2/sum(x$model)}
results[i,"om.model"] <- omega.model
observed.cor <- cor(x$observed)
ev <- eigen(observed.cor)$values
f1 <- fa(observed.cor)$loadings
om.fa <- sum(f1)^2/sum(observed.cor)
e.f1 <- sum(f1^2)/nvar
    sem.model <- omegaSem(x$fload[,1:nfact],sl=TRUE,nfactors=nfact)  #this is the model based upon the true values
    if(sem) {stop('The option to use the sem package  has been replaced with calls to lavaan')
    #if(!requireNamespace('sem')) {stop("You must have the sem package installed to use omegaSem")} else {sem.om <- try(sem(model=sem.model,S=observed.cor, N=n))} 
 
  #omega.cfa <- omegaFromSem(observed.cor,sem.om,flip=flip)
  #  if(omega.cfa$omega >1) omega.cfa$omega <- NA
  #results[i,"omegaCFA"] <- omega.cfa$omega
  } else {omega.cfa <- NULL}
results[i,"n"] <- n
if(n > 0) {
   if (sem) {om <- try(omegaSem(x$observed,om.fact,flip=flip,plot=FALSE,option=option))} else {
om <- try(omega(x$observed,om.fact,flip=flip,plot=FALSE,option=option))}
          ic <- suppressWarnings(ICLUST(x$observed,1,plot=FALSE))} else {
           if (sem) {om <- try(omegaSem(x$model,om.fact,flip=flip,plot=FALSE,option=option))} else {om <- try(omega(x$model,om.fact,flip=flip,plot=FALSE,option=option))
           if(inherits(om,"try-error")) {message("Error in sem. iteration = ",i)
            om <- NA
            next}}
          ic <- suppressWarnings(ICLUST(x$model,1,plot=FALSE))}
#results
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
  
  
#not public
#simulate the difference between  two groups  
  sim.groups <- function(n=1000,r=.5,d=.5,nvar = 2,bg= -1) {
model <- matrix(r,nvar,nvar)
 diag(model) <- 1
 eX <- eigen(model)
 observed <- matrix(rnorm(nvar * n),n)
 observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed))
 n1 <- n/2
if(bg < 1 ) {  mu1 <- c(0,d)} else {mu1 <- 0 }
 group1 <- t(t( observed[1:n1,]) + mu1  )
if(bg < 1) {mu2 <- c(d,0)} else {mu2 <- d}
 group2 <- t( t(observed[(n1+1):n,]) + mu2)
 data <- data.frame(grp=1,vars=group1)
 data2 <- data.frame(grp=2,vars=group2)
 data <- rbind(data,data2)
 return(data)
 }
 
 
 "sim.general.theta" <- 
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