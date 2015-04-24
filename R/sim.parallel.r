"sim.parallel" <-
function(ntrials=10,nvar = c(12,24,36,48),nfact = c(1,2,3,4,6),
n = c(200,400)) {
nvariables = nvar
factors = nfact
subjects = n 

result <- matrix(NaN,ncol=7,nrow=ntrials*length(nvariables) * length(subjects) * length(factors))
k <- 1
for (nfact in factors) {
  for (nvar in nvariables) {
    for (nsub in subjects) {
   for (trials in 1:ntrials) {
   x <- sim.minor(nvar=nvar,nfact=nfact,n=nsub)$observed
   fp <- fa.parallel(x)
   fps <- fa.parallel(x,SMC=TRUE) 
   result[k,1] <- nfact
   result[k,2] <- nvar
   result[k,3] <- trials
   result[k,4] <- fp$nfact
   result[k,5] <- fps$nfact
   result[k,6] <- fp$ncomp
   result[k,7] <- nsub
   k <- k + 1 
  } #trials
 } #subjects
 }#variables
}#factors
colnames(result) <- c("factors","nvar","trials","nfact","smc.fact","ncomp","nsub")
return(result)
}

"sim.correlation" <- function(R,n=1000,data=FALSE) {
     eX <- eigen(R)
     nvar <- ncol(R) 
     observed <- matrix(rnorm(nvar * n),n,nvar)
     observed <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed))
    colnames(observed) <- colnames(R)
    if(data) {result <- observed} else {
  	 result <- cor(observed)}
  	 return(result)}
