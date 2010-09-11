"omegaSem" <-
function(m,nfactors=3,fm="minres",key=NULL,flip=TRUE, digits=2,title="Omega",sl=TRUE,labels=NULL, plot=TRUE,n.obs=NA,rotate="oblimin",Phi = NULL,option="equal",...) {
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      #if Phi is not null, this implies that we have been given a factor matrix  -- added May 30, 2010
if(!require(sem)) stop("You must have the sem package installed to use omegaSem")
if (!sl) {warning("OmegaSem only works for Bifactor models, sl set to TRUE ")
    sl <- TRUE}

    cl <- match.call() 
om <- omega(m,nfactors,fm,key,flip, digits,title,sl,labels, plot,n.obs,rotate,Phi,option,...)
      #m is a correlation matrix, or if not, the correlation matrix is found
      #nfactors is the number of factors to extract
      #key allows items to be reversed scored  if desired
      #if Phi is not null, this implies that we have been given a factor matrix  -- added May 30, 2010
      
sem.model <- om$model
if(is.na(n.obs)) {warning("n.obs must be specified to use omegaSem.  n.obs arbitrarily set to 500")
               n.obs <- 500}

 if(dim(m)[1] != dim(m)[2]) {
                            n.obs <- dim(m)[1]
                            m <- cor(m,use="pairwise")} else {
                            m <- cov2cor(as.matrix(m))    #make sure it is a correlation matrix not a covariance or data matrix (if we change this, we will need to change the calculation for omega later)
                           }
      nvar <- dim(m)[2]
if(is.na(n.obs)) {message("Number of observations not specified. Arbitrarily set to 500")
                 n.obs <- 500
                 }
                            
      if(is.null(colnames(m))) {  rownames(m) <- colnames(m) <- paste("V",1:nvar,sep="") }
       m.names <- colnames(m)
if(!require(sem)) {stop("You must have the sem package installed to use omegaSem")
 } else {sem.om <- sem(sem.model,m, n.obs) }
omega.efa <- omegaFromSem(m,sem.om,flip=flip)
results <- list(omegaSem=om,omega.efa=omega.efa,sem=sem.om,Call=cl)
class(results) <- c("psych", "omegaSem")
return(results)
}


"omegaFromSem" <- 
function(m,s,flip=TRUE) {
# m is the correlation matrix
# s is the sem solution
nvar <- dim(m)[1]
n.fact <-dim(s$A)[2]-nvar
rn <- rownames(m)
g <- s$A[1:nvar,nvar+n.fact]


cfa.loads <- cbind(g,s$A[1:nvar,(nvar+1):(nvar+n.fact-1)])  #this puts g last

if(flip) { 
	flipper <- rep(1,nvar)
	flipper[g < 0] <- -1
	signkey <- strtrim(flipper,1)
            		 signkey[signkey=="1"] <- ""
            		 rn <- paste(rn,signkey,sep="")
	flipper <- diag(flipper)
	m <- flipper %*% m %*% t(flipper) 
	g <- flipper %*% g
	cfa.loads <- flipper %*% cfa.loads}


w <- solve(m,cfa.loads)
rownames(cfa.loads)  <- rn
rownames(w)  <- rn

gR2 <- diag(t(w) %*% cfa.loads)

om.sem <- sum(g)^2/sum(m)
class(cfa.loads) <- "loadings"
results <- list(omega=om.sem,cfa.loads=cfa.loads,gR2=gR2)
return(results)
}


