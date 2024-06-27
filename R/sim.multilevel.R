sim.multilevel <- function(nvar=9,ngroups=4,ncases=16,rwg,rbg,eta) {
 e.wg <- eigen(rwg)
 v.wg <- pmax(e.wg$values,0)
 etabg <- sqrt(1-eta^2)
 e.bg <- eigen(rbg)
 v.bg <- pmax(e.bg$values,0)
 wg<- matrix(rnorm(nvar*ncases),ncases)
 wg <- scale(wg)
 wg <- t(e.wg$vectors %*% sqrt(diag(v.wg)) %*% t(wg))
 bg <- matrix(rnorm(nvar*ngroups),ngroups)
 bg <- scale(bg)
 bg <- e.bg$vectors %*% sqrt(diag(v.bg)) %*% t(bg)
 bg <- matrix(rep(bg, (ncases/ngroups)),nrow=ncases,byrow=TRUE)
 gr <- rep((1:ngroups),(ncases/ngroups))
 XY <- wg %*% diag(eta^2)   +  bg %*% diag(etabg^2) 
 XY <- cbind(gr,XY)
 colnames(XY) <- c("Group",paste("V",1:nvar,sep=""))
 result <- list(wg=wg,bg=bg,xy=XY)
 }
 
 
 
 #Created January 28, 2017
#meant to simulate a number of within subject multilevel models


#Created January 28, 2017
#meant to simulate a number of within subject multilevel models


"sim.multi" <-
 function(n.obs=4,nvar = 2,nfact=2, ntrials=96,days=16,mu=0,sigma=1,fact=NULL,loading=.9,phi=0,phi.i = NULL,beta.i=0,mu.i=0,sigma.i = 1,sin.i=0,cos.i=0,AR1=0,f.i=NULL,plot=TRUE) {
if(missing(n.obs)) n.obs=4
X <- list()
Xjk <- matrix(NA,ncol=nvar +nfact,nrow=ntrials)
if(missing(mu) ) mu <- rep(0,nvar) 
if(missing(sigma)) sigma <- rep(1,nvar)
#if(missing(phi.i)) {phi.i <-  phi} 
if(missing(beta.i)) { beta.i <- matrix(0,ncol=nvar,nrow=n.obs) } else {
     if(length(beta.i) < n.obs) {beta.i <- matrix(beta.i,ncol=nvar,nrow=n.obs,byrow=TRUE) } }
if(missing(mu.i)) mu.i <- matrix(0,ncol=nvar,nrow=n.obs)
if(missing(sigma.i)) {sigma.i <- matrix(1,ncol=nvar,nrow=n.obs)} else {
        if(length(sigma.i) < n.obs) {sigma.i <- matrix(sigma.i,ncol =nvar,nrow=n.obs,byrow=TRUE)}
       }
if(missing(sin.i)) {sin.i <- matrix(0,ncol=nvar,nrow=n.obs)} else {
  if (length(sin.i) < n.obs) {sin.i <- matrix(sin.i,ncol=nvar,nrow=n.obs,byrow=
  TRUE) }
  }
if(missing(cos.i)) {cos.i <- matrix(0,ncol=nvar,nrow=n.obs) } else {
  if (length(cos.i) < n.obs) {cos.i <- matrix(cos.i,ncol=nvar,nrow=n.obs,byrow=
  TRUE) }
  }
if(missing(AR1)) {AR1 <- matrix(0,ncol=nvar,nrow=n.obs) }  else {
  if (length(AR1) < n.obs) {AR1 <- matrix(AR1,ncol=nfact,nrow=n.obs,byrow=
  TRUE) }
  }
  
  if(is.null(phi)) {phi <-diag(1,nfact) } else {phi <- matrix(phi,ncol=nfact,nrow=nfact)
              diag(phi) <- 1}
              
  
    
    
  if(!is.null(phi.i)) {if(length(phi.i) < n.obs) {phi.i <- rep(phi.i,n.obs/length(phi.i))} }  
  if(nfact > 1) {

  if(is.null(fact)) {     #these are the group level factor loadings
   fact <- matrix(0,nrow=nvar,ncol=nfact)
   # fact[ ] =((( col(fact)+ row(fact)) %% nfact ))  * loading
  # fact[((round(row(fact)/nvar))+1) == col(fact)] <- loading  #just works for two factors!
  for(j in 1:nfact) {
    #fact[((j-1)*nvar/nfact +1):j*nvar/nfact,j] <- loading
    fact[((j-1) * nvar/nfact +1):(j*nvar/nfact),j] <- loading
   }  
    
    fact<-  (fact %*% phi )}} else { fact <- matrix(loading,ncol=nfact,nrow=nvar) }
 if(is.null(f.i)) { f.i <- list()
 for (i in 1:n.obs) {
   f.i[[i]] <- fact
     }
     }


trials.day <- ntrials/days
hours <- 24/trials.day
time <- seq(hours,days * trials.day*hours,hours)
t.radian <- time * pi /12


for (i in 1:n.obs)  {
  xij <- rnorm((nvar + nfact),mu,sigma)   #between subjects 
for(j in 1:nfact) {   #first generate the factor scores that have a within subject model
  error <- rnorm(ntrials,mu.i[i,j],sigma.i[i,j])
  lagerror <- c(0, error[1:(ntrials-1)])
  Xjk[,j] <- xij[j] +  mu[j] +beta.i[i,j] *time/ntrials + sin(t.radian)*sin.i[i,j] + cos(t.radian)*cos.i[i,j] +  error + AR1[i,j] * lagerror
  
  }
  
  #factor scores are the first nfact elemements of Xjk
  #now, generate item scores
  

  if(is.null(phi.i)) {phi.ind <- diag(1,nfact) } else {phi.ind <- matrix(phi.i[i],nfact,nfact)
               diag(phi.ind) <- 1}
               
  Xjk[,1:nfact] <- Xjk [,1:nfact] %*% phi.ind     #these are the factor scores for the ith subject



for(k in 1:nvar) {
   uniq <- sqrt(1 - sum(f.i[[i]][k,]^2))   #the uniqueness is 1-h2
   uniq.err <-  rnorm(ntrials,0,uniq)
   score <- 0
      for (j in 1:nfact) {
         score <- score +  Xjk[,j] * f.i[[i]][k,j]    #these are orthogonal factor scores  * loadings    -- can add phi into this 
        }
    Xjk[,nfact  + k ] <- score + uniq.err  #these are the variable scores
    }
X[[i]] <- Xjk    #save them for this subject
}

#now, lets summarize what we have done
DV <- unlist(X)  

#This organizes the data with a separate column for each variable for analysis by statsBy
dv.a <- array(DV,dim=c(ntrials,nvar+ nfact,n.obs))
dv.m <-NULL
for(i in 1:(nvar+nfact)) { dv.m <- cbind(dv.m,as.vector(dv.a[,i,])) }

dv.df <- data.frame(dv.m,time = rep(time,n.obs),id=rep(1:n.obs,each=ntrials))
colnames(dv.df)[1:(nfact+nvar)] <- c(paste0("F",1:nfact),paste0("V",1:nvar))

#
#However, if we want to plot it, we need to rearrange a bit more
if(plot) {
IV <- NULL #to get around the problem that RCMD check thinks this is global
vars <-  c(paste0("F",1:nfact),paste0("V",1:nvar))

select <- rep(NA,nvar * ntrials * n.obs) #we use select to pick out the variables, but leave the factors
kount <- ntrials*nfact
for(i in 1:n.obs) {
    select[(1:(ntrials*nvar)+(i-1) * ntrials *nvar)]  <- kount + 1:(ntrials*nvar)
    kount <- kount + ntrials *( nvar+nfact)
    }

vars <- paste0("V",1:nvar)
X.df <- data.frame(DV = DV[select], time=rep(time,(n.obs*(nvar))),id = rep(1:n.obs,each=(ntrials*(nvar))),IV =  rep(rep(vars,each=ntrials),n.obs) )

plot1<- xyplot(DV ~ time | id, group=IV, data=X.df, type = "b",as.table=TRUE,strip=strip.custom(strip.names=TRUE,strip.levels=TRUE),col=c("blue","red","black","grey"))
print(plot1) }

invisible(dv.df)   #this returns the scores if we want to do further processing
}


 