#Added June 20, 2018 to try to do Neils Waller's direct Schmid Leiman

Procrustes <-function(L, Target=NULL){#Adapted from Niels Waller (2017)
if(is.null(Target)) Target <- factor2cluster(L) 
tM1M2 <- t(Target) %*% L
 svdtM1M2 <- svd(tM1M2) 
 P <- svdtM1M2$u
 Q <- svdtM1M2$v 
T <- Q%*%t(P)
## Orthogonally rotate L to Target
 return(list(loadings = L %*%T,rotation = T))
 }
 
 #allowing to specify a number of rotations
 oblique.rotations <- function(rotate="oblimin",loadings,...){
   if (rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax" | rotate =="geominQ"  | rotate =="bentlerQ"  |rotate == "targetQ"  ) {
     				if (!requireNamespace('GPArotation')) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
     				    Phi <- NULL} else { 
     				      
     				             ob <- try(do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...)))
     				               if(class(ob)== as.character("try-error"))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
     				               ob <- Promax(loadings)}
     				                 
     				loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}
  		                             }
 return(list(loadings=loadings,Phi=Phi))
 }
 
 
 #direct Schmid Leiman adapted from Waller (2017)
 direct.sl <- function(x,nfactors=3,rotate="oblimin",cut=.3,plot=TRUE){
 nvar <- ncol(x)
 f <- fa(x,nfactors=nfactors,rotate ='none')   #unrotated solution
 #construct the target from the rotated solution
 f.obl <- oblique.rotations(rotate=rotate,loadings = f$loadings)$loadings
 targ <- factor2cluster(f.obl,cut=cut)
 #Waller adjustments to target and factor model
 targ <- cbind(g=rep(1,nvar),targ)
 f0 <- cbind(rep(0,nvar),f$loadings)
 direct <- Procrustes(f0,targ)$loadings   #The Waller Procrustes solution
 colnames(direct) <- colnames(targ)  #put some labels in
 fa.diagram(direct,g=TRUE,simple=FALSE,cut=cut)
 return(direct)
 }
