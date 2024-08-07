
#modified 7/25/21 to report factor scores so that we can use biplot on the exensions.
#modified again 12/03/23 to report scores as scores rather than with all the extra baggage
"fa.extension" <-
  function(Roe,fo,correct=TRUE) {
 cl <- match.call()
 omega <-FALSE

 if(!is.null(class(fo)[2])) {if(inherits(fo,"fa")) {
          if(!is.null(fo$Phi)) {Phi <- fo$Phi} else {Phi <- NULL}
          
       fl <- fo$loadings 
       fs <- fo$Structure
       
     } else {if (inherits(fo,"omega")) {    #switched to inherits  December 20, 2019
         omega <- TRUE
         w <- fo$stats$weights
         fl <- fo$schmid$sl
         Phi <- NULL
         fl <- fl[,1:(dim(fl)[2]-4)]    #just the loadings
         nfactors <- dim(fl)[2]
         fe <- t(t(w) %*% Roe)
         foblique <- fo$schmid$oblique
         feoblique <- t( Roe) %*% foblique %*% (solve(t(foblique)%*% (foblique)))
         feoblique <- feoblique %*% solve(fo$schmid$phi)  
       } 
    }
    }
 #Roe is Horn's Re    R1 is Phi  Pc is pattern of original = fl
 # Pe = Re Pc  solve (Pc'Pc) solve Phi
 if(!omega) fe <- t( Roe) %*% fl %*% (solve(t(fl)%*% (fl)))  #should we include Phi?
  if(!is.null(Phi)) fe <- fe %*% solve(Phi)  #horn equation 26
 
 if(!correct) {#the Gorsuch case -- not actually-- read Gorsuch again

    # d <- diag(t(fl) %*% fo$weight)  #this is probably  wrong 
     
     d <- sqrt(diag(t(fl) %*% fo$weight))   #a correction of sorts for reliability
     fe <- (fe * d)
 }
 colnames(fe) <- colnames(fl)
 rownames(fe) <- colnames(Roe)
if(!is.null(Phi)) {resid <- Roe - fl %*% Phi %*% t(fe)} else {resid <- Roe - fl  %*% t(fe)}   #fixed to actually give residual  (1/30/18)
 result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)
 if(!omega)  {result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)} else {result <- list(loadings = fe,oblique= feoblique,Phi=Phi,resid=resid,Call=cl)}
 class(result) <- c("psych","extension")
 return(result)
}
#written April 5, 2011
#revised August 15, 2011 to avoid using the weights matrix except in the omega case

#created December 8, 2012 to allow for extension and goodness of fits of total model
#modified 31/5/14 to allow for omega extension as well 
#modified 04-09/16 to pass the Structure matrix as well
#Added the cors and correct parameters to pass to fa 1/3/21


"fa.extend" <- 
function(r,nfactors=1,ov=NULL,ev=NULL,n.obs = NA, np.obs=NULL,correct=TRUE,rotate="oblimin",SMC=TRUE,warnings=TRUE, fm="minres",alpha=.1, omega=FALSE,cor="cor",use="pairwise",cor.correct=.5,weight=NULL,missing=FALSE,smooth=TRUE, ...) {
 cl <- match.call()
 if(is.numeric(ev)) ev  <- colnames(r)[ev]    #in case we are selecting variables 
 if(is.numeric(ov)) ov  <- colnames(r)[ov]
  nv <- c(ov,ev)
 #if(nrow(r) > ncol(r)){  #the case of a data matrix
 if(!isCorrelation(r)){
     #first find the correlations
      n.obs <- nrow(r)
         np.obs.r <- pairwiseCount(r)[nv,nv]
         np.obs <- np.obs.r[ov,ov]
         data <- r  #if we want to find factor scores
  # r <- cor(r,use='pairwise') 
       switch(cor, 
       cor = {r <- cor(r,use=use) },   #does not include the weight option from fa 
       cov = {r <- cov(r,use=use) 
              covar <- TRUE},
       wtd = { r <- cor.wt(r,w=weight)$r},
       spearman = {r <- cor(r,use=use,method="spearman")},
       kendall = {r <- cor(r,use=use,method="kendall")},
       tet = {r <- tetrachoric(r,correct=cor.correct,weight=weight)$rho},
       poly = {r <- polychoric(r,correct=cor.correct,weight=weight)$rho},
       tetrachoric = {r <- tetrachoric(r,correct=cor.correct,weight=weight)$rho},
       polychoric = {r <- polychoric(r,correct=cor.correct,weight=weight)$rho},
       mixed = {r <- mixedCor(r,use=use,correct=cor.correct)$rho}
       )
       
 if(omega) {fo <- omega(r[ov,ov],nfactors=nfactors,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,...)} else {
       fo <- fa(r[ov,ov],nfactors=nfactors,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,cor=cor,alpha=alpha, missing=missing,impute="mean", smooth=smooth,...)}
     
           
    } else {  #the case of a correlation matrix  
       data <- NULL       
       R <- r[ov,ov]
       np.obs.r <- np.obs
      if(omega) {fo <- omega(R,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,cor=cor,alpha=alpha,np.obs=np.obs[ov,ov],...)} else { 
      fo <- fa(R,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,cor=cor, correct=correct,alpha=alpha,np.obs=np.obs[ov,ov],missing=missing,impute="mean",smooth=smooth,...)}
     }
Roe <- r[ov,ev,drop=FALSE]
fe <- fa.extension(Roe,fo,correct=correct)
if(omega) fo$loadings <- fo$schmid$sl[,1:(ncol(fo$schmid$sl)-3)]
       
foe <- rbind(fo$loadings,fe$loadings)

if(omega) oblique <- rbind(fo$schmid$oblique,fe$oblique)

if(is.na(n.obs) && !is.null(np.obs)) n.obs <- max(as.vector(np.obs))
result <- factor.stats(r[nv,nv],foe,fo$Phi,n.obs,np.obs.r,alpha=alpha,smooth=smooth,coarse=FALSE)  #this is the second call to factor.stats
if(omega) result$schmid$sl <- foe
    result$rotation <- rotate
    result$loadings <- foe
    if(nfactors > 1) {if(is.null(fo$Phi)) {h2 <- rowSums(foe^2)} else {h2 <- diag(foe %*% fo$Phi %*% t(foe)) }} else {h2 <-foe^2}
    result$communality <- h2
    result$fm <- fm  #remember what kind of analysis we did
    
    result$fo=fo
    if(!is.null(data)) result$scores <- factor.scores(data[,ov],fo,missing=missing,impute="mean")$scores
    if(omega) {result$schmid$sl <- foe
              result$schmid$gloading <- fo$schmid$gloading
              result$schmid$oblique <- oblique
            }
   if(is.null(fo$Phi)) {result$Structure <- foe } else { result$Structure <- foe %*% fo$Phi}
    result$fe=fe
    result$resid=fe$resid
    result$Phi=fo$Phi
    result$fn="fa"
    result$Call=cl

class(result) <- c("psych","extend")
return(result)
}



