

"fa.extension" <-
  function(Roe,fo,correct=TRUE) {
 cl <- match.call()
 omega <-FALSE

 if(!is.null(class(fo)[2])) {if(class(fo)[2]=="fa") {
          if(!is.null(fo$Phi)) {Phi <- fo$Phi} else {Phi <- NULL}
          
       fl <- fo$loadings 
       
     } else {if (class(fo)[2] == "omega") {
         omega <- TRUE
         w <- fo$stats$weights
         fl <- fo$schmid$sl
         Phi <- NULL
         fl <- fl[,1:(dim(fl)[2]-3)]
         nfactors <- dim(fl)[2]
         fe <- t(t(w) %*% Roe)
         foblique <- fo$schmid$oblique
         feoblique <- t( Roe) %*% foblique %*% (solve(t(foblique)%*% (foblique)))
         feoblique <- feoblique %*% solve(fo$schmid$phi)  
       } 
    }
    }
 
 if(!omega) fe <- t( Roe) %*% fl %*% (solve(t(fl)%*% (fl))) 
  if(!is.null(Phi)) fe <- fe %*% solve(Phi)
 
 if(!correct) {#the Gorsuch case
     d <-diag(t(fl) %*% fo$weight)
     fe <- (fe * d)
 }
 colnames(fe) <- colnames(fl)
if(!is.null(Phi)) {resid <- Roe - fl %*% Phi %*% t(fe)} else {resid <- fl  %*% t(fe)}
 result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)
 if(!omega)  {result <- list(loadings = fe,Phi=Phi,resid=resid,Call=cl)} else {result <- list(loadings = fe,oblique= feoblique,Phi=Phi,resid=resid,Call=cl)}
 class(result) <- c("psych","extension")
 return(result)
}
#written April 5, 2011
#revised August 15, 2011 to avoid using the weights matrix except in the omega case

#created December 8, 2012 to allow for extension and goodness of fits of total model
#modified 31/5/14 to allow for omega extension as well 
"fa.extend" <- 
function(r,nfactors=1,ov=NULL,ev=NULL,n.obs = NA, np.obs=NULL,correct=TRUE,rotate="oblimin",SMC=TRUE,warnings=TRUE, fm="minres",alpha=.1, omega=FALSE, ...) {
 cl <- match.call()
  nv <- c(ov,ev)
 if(nrow(r) > ncol(r)){  #the case of a data matrix
 if(omega) {fo <- omega(r[,ov],nfactors=nfactors,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,...)} else {
       fo <- fa(r[,ov],nfactors=nfactors,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,...)}
         n.obs <- nrow(r)
         np.obs.r <- count.pairwise(r)[nv,nv]
         np.obs <- np.obs.r[ov,ov]
        r <- cor(r,use='pairwise')  
    } else {  #the case of a correlation matrix         
       R <- r[ov,ov]
       np.obs.r <- np.obs
      if(omega) {fo <- omega(R,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,np.obs=np.obs[ov,ov],...)} else { 
      fo <- fa(R,nfactors=nfactors,n.obs=n.obs,rotate=rotate,SMC=SMC,warnings=warnings,fm=fm,alpha=alpha,np.obs=np.obs[ov,ov],...)}
     }
Roe <- r[ov,ev]
fe <- fa.extension(Roe,fo,correct=correct)
if(omega) fo$loadings <- fo$schmid$sl[,1:(ncol(fo$schmid$sl)-3)]
       
foe <- rbind(fo$loadings,fe$loadings)

if(omega) oblique <- rbind(fo$schmid$oblique,fe$oblique)

if(is.na(n.obs) && !is.null(np.obs)) n.obs <- max(as.vector(np.obs))
result <- factor.stats(r[nv,nv],foe,fo$Phi,n.obs,np.obs.r,alpha=alpha)
if(omega) result$schmid$sl <- foe
    result$rotation <- rotate
    result$loadings <- foe
    if(nfactors > 1) {if(is.null(fo$Phi)) {h2 <- rowSums(foe^2)} else {h2 <- diag(foe %*% fo$Phi %*% t(foe)) }} else {h2 <-foe^2}
    result$communality <- h2
    result$fm <- fm  #remember what kind of analysis we did
    
    result$fo=fo
    if(omega) {result$schmid$sl <- foe
              result$schmid$gloading <- fo$schmid$gloading
              result$schmid$oblique <- oblique
            }
    result$fe=fe
    result$Phi=fo$Phi
    result$fn="fa"
    result$Call=cl
class(result) <- c("psych","extend")
return(result)
}

 
 
