#corrected estimate of communality, May 21, 2007
#removed "x" from factanal call, June 14, 2007
#added ability to do 2 factors by treating them with equal loadings Jan 2008
#added use of simplimax rotation June 2008
#corrected the sign of group factors to match the original factors
#modified November, 2014 to allow for covariances in addition to correlations.
#also cleaned up the code to take advantage of switch
"schmid" <-
function (model, nfactors = 3, fm = "minres",  digits=2,rotate="oblimin",n.obs=NA,option="equal",Phi=NULL,covar=FALSE,two.ok=FALSE,...) 
{
 cl <- match.call()
#if Phi is not Null, then we have been given a factor matrix, otherwise
#model is a correlation matrix, or if not, the correlation matrix is found
#nfactors is the number of factors to extract
if(!requireNamespace('GPArotation')) {stop("I am sorry, you need to have the  GPArotation package installed")}
#there are two types of input
# 1: from a factor analysis and or rotation function with a matrix of loadings and a Phi matrix
# 2: raw data or a correlation/covariance matrix  (from e.g, omega)
# 3 Input is the output of a factor analysis 
  if(is.list(model)) {if((class(model)[1] == "psych") && (class(model)[2] == "fa")) {Phi <- model$Phi
                model <- model$loadings} else {stop("input is a list, but is not from of class 'fa' ")}}
  if(is.null(Phi)) {  #the normal case
      normal.case <- TRUE
      nvar <-dim(model)[2]
      if(dim(model)[1] != dim(model)[2]) {n.obs <- dim(model)[1]
                                        if(covar) { model <- cov(model,use="pairwise")} else {model <- cor(model,use="pairwise") }
                                         }
                                          
     if (fm =="pc") {
        fact <- principal(model, nfactors,n.obs=n.obs,covar=TRUE,...)
        fm <- 'minres'    #because we want to use factor analysis for the higher level factors
        message("The higher order factor is found using minres -- see the notes")
    } else {if ((fm == "pa") |(fm =="minres") | (fm =="wls")  |(fm =="minres") |(fm =="ml")|(fm =="mle")  |(fm =="gls") |(fm =="minchi") |(fm =="minrank")) {fact <- fa(model, nfactors,n.obs=n.obs,rotate="varimax",fm=fm,covar=covar) } else {
     
        stop("The method of factor extraction you specified is not available")
        
    }}
     orth.load <- loadings(fact)
    } else {model <- as.matrix(model)
            Phi <- as.matrix(Phi)
            fact <- model %*% Phi  #find the orthogonal (structure) matrix from the oblique pattern and the Phi matrix
            orth.load <- fact   #but this is not the correct factor  solution
            ev <- eigen(Phi)
            orth.load <- model %*% ev$vector %*% sqrt(diag(ev$values))
            colnames(orth.load) <- colnames(Phi)
            nfactors <- dim(fact)[2]
            normal.case <-FALSE}
   
    colnames(orth.load)  <- paste("F",1:nfactors,sep="")
    if(nfactors == 1) { message("Omega_h for 1 factor is not meaningful, just omega_t")
                        obminfact <-list(loadings= orth.load)
                       factr <- 1
                       
           } else {  #the normal case is nfactors > 2
switch(rotate,
    simplimax = {obminfact <- GPArotation::simplimax(orth.load)},
 #    promax  =    {obminfact  <- Promax(orth.load)
#      								 rotmat <- obminfact$rotmat
#                    						Phi <- obminfact$Phi
#            							 },
       promax =   {#to match revised fa call 
                pro <- kaiser(orth.load,rotate="Promax",...)   #calling promax will now do the Kaiser normalization before doing Promax rotation
     			 obminfact <- pro
     			  rot.mat <- pro$rotmat
     			  Phi <- pro$Phi 
     			},	      							 
    Promax = {obminfact  <- Promax(orth.load)
     								 rotmat <- obminfact$rotmat
                   						Phi <- obminfact$Phi
           							 },
    TargetQ = {obminfact <- do.call(rotate,list(orth.load,...)) 
           							       loadings <- obminfact$loadings
     			                            Phi <- obminfact$Phi
     			                       },
     			                       
     cluster  = {obminfact <- varimax(orth.load)            			
								obminfact <- target.rot(obminfact,...)
     			              	loadings <- obminfact$loadings
     			                Phi <- obminfact$Phi
     			                 },
     target = {obminfact <- varimax(orth.load)            			
								obminfact <- target.rot(obminfact,...)
     			              	loadings <- obminfact$loadings
     			                Phi <- obminfact$Phi
     			                 },
     oblimin = 	{ obminfact <- try(GPArotation::oblimin(orth.load))
           							        if(inherits(obminfact,as.character("try-error"))) {obminfact <- Promax(orth.load)   #special case for examples with exactly 2 orthogonal factors
           							        message("\nThe oblimin solution failed, Promax used instead.\n")                   #perhaps no longer necessary with patch to GPForth and GPFoblq in GPArotation
           							        rotmat <- obminfact$rotmat
                   						    Phi <- obminfact$Phi}},
       geominQ = 	{ obminfact <- try(GPArotation::geominQ(orth.load))
           							        if(inherits(obminfact,as.character("try-error"))) {obminfact <- Promax(orth.load)   #special case for examples with exactly 2 orthogonal factors
           							        message("\nThe geominQ solution failed, Promax used instead.\n")                   #perhaps no longer necessary with patch to GPForth and GPFoblq in GPArotation
           							        rotmat <- obminfact$rotmat
                   						    Phi <- obminfact$Phi}},
                   						    
       bentlerQ = 	{ obminfact <- try(GPArotation::bentlerQ(orth.load))
           							        if(inherits(obminfact,as.character("try-error"))) {obminfact <- Promax(orth.load)   #special case for examples with exactly 2 orthogonal factors
           							        message("\nThe bentlerQ solution failed, Promax used instead.\n")                   #perhaps no longer necessary with patch to GPForth and GPFoblq in GPArotation
           							        rotmat <- obminfact$rotmat
                   						    Phi <- obminfact$Phi}},            						    
        targetQ =   { obminfact <- try(GPArotation::targetQ(orth.load,...))
           							        if(inherits(obminfact,as.character("try-error"))) {obminfact <- Promax(orth.load)   #special case for examples with exactly 2 orthogonal factors
           							        message("\nThe targetQ solution failed, Promax used instead.\n")                   #perhaps no longer necessary with patch to GPForth and GPFoblq in GPArotation
           							        rotmat <- obminfact$rotmat
                   						    Phi <- obminfact$Phi}},          						    
                   						      
         biquartimin =    {obminfact <- biquartimin(orth.load,...)
                    loadings <- obminfact$loadings
     				 Phi <- obminfact$Phi
     				 rot.mat <- t(solve(obminfact$Th))}
                   						    
               						    
     )
                        
     			                       
     			                       
     			                       
#     			       
#           						    
#      if (rotate == "simplimax") {obminfact <- simplimax(orth.load)} else {
#      if((rotate == "promax") | (rotate == "Promax")  )    {obminfact  <- Promax(orth.load)
#     								 rotmat <- obminfact$rotmat
#                   						Phi <- obminfact$Phi
#           							 } else {
#           							        if(rotate=="TargetQ") {obminfact <- do.call(rotate,list(orth.load,...)) 
#           							       loadings <- obminfact$loadings
#     			                            Phi <- obminfact$Phi
#     			                       } else {
#           							     
#           							 if ((rotate == "cluster") | (rotate == "target")) {obminfact <- varimax(orth.load)            			
#								obminfact <- target.rot(obminfact,...)
#     			              	loadings <- obminfact$loadings
#     			                Phi <- obminfact$Phi
#     			                 } else {
#           							  obminfact <- try(oblimin(orth.load))
#           							        if(class(obminfact)== as.character("try-error")) {obminfact <- Promax(orth.load)   #special case for examples with exactly 2 orthogonal factors
#           							        message("\nThe oblimin solution failed, Promax used instead.\n")                   #perhaps no longer necessary with patch to GPForth and GPFoblq in GPArotation
#           							        rotmat <- obminfact$rotmat
#                   						    Phi <- obminfact$Phi} }} }
#                   						 }
          		}  
    if(nfactors > 1) rownames(obminfact$loadings) <- attr(model,"dimnames")[[1]]
    
    if(!normal.case) { fload <- model
                        factr <- Phi
                        model <- fload %*% Phi %*% t(fload)
                        diag(model) <- 1} else {
                    	fload <- obminfact$loadings
                                #factr <- t(obminfact$Th) %*% (obminfact$Th)
                    	factr <- obminfact$Phi}
   
   if (nfactors ==1) {gload <- c(1)
              warning("Omega_h and Omega_asymptotic are not meaningful with one factor") } else { colnames(factr) <- rownames(factr) <- paste("F",1:nfactors,sep="")  #make it a vector
   if (nfactors>2) {
   
     
     gfactor <- fa(factr,fm=fm)   #The first factor of the factor intercorrelation matrix
                  #added fm=fm  March 5, 2011
    
       gload <- loadings(gfactor) } else {gload<- c(NA,NA)   #consider the case of two factors 
            if(option=="equal") {
      			 gload[1] <- sqrt(abs(factr[1,2]))
      			 gload[2] <- sign(factr[1,2])*sqrt(abs(factr[1,2])) 
      			if(!two.ok) message("\nThree factors are required for identification -- general factor loadings set to be equal. \nProceed with caution. \nThink about redoing the analysis with alternative values of the 'option' setting.\n")} else { if(option=="first") {
      			 gload[1] <- 1
      			# gload[2] <- abs(factr[1,2])
      			gload[2] <-  (factr[1,2])
      			 message("\nThree factors are required for identification -- general factor loading set to be 1 for group factor 1. \nProceed with caution. \nThink about redoing the analysis with alternative values of the 'option' setting.\n")} else { gload[2] <- 1
      			# gload[1] <- abs(factr[1,2]) 
      			 gload[1] <- (factr[1,2]) 
      			 message("\nThree factors are required for identification -- general factor loadings are set to be 1 for group factor 2.\nProceed with caution. \nThink about redoing the analysis with alternative values of the 'option' setting.\n")} }
      			 
       
              }  
    }
    gprimaryload <- fload %*% gload
    colnames(gprimaryload) <- "g"
    h2 <- diag(orth.load %*% t(orth.load)) 
   # u2 <- 1 - h2    
    u2 <- diag(model) - h2                     
    uniq <- diag(model)- fload^2
    guniq <- diag(model) - gprimaryload^2
    #Ig <- matrix(0, ncol = nfactors, nrow = nfactors)
    #diag(Ig) <- gload
    Ig <- diag(drop(gload))   #3/5/11
    primeload <- fload %*% Ig
    g.percent <- gprimaryload^2/h2
    colnames(g.percent) <- "p2"
    uniq2 <- diag(model) - uniq - primeload^2
    uniq2[uniq2<0] <- 0
    sm <-  sign(fload) * sqrt(uniq2)  #added June 1, 2010 to correctly identify sign of group factors
    F <- cbind(gprimaryload, sm)  #the factor pattern matrix
    #the next two lines are actually not needed because the factors are orthogonal
    Structure  <- t( Pinv(t(F) %*% F) %*% t(F) %*% orth.load %*% t(orth.load))
    Phi.S <- t(Structure) %*% F %*% Pinv(t(F) %*% F)  #this is the pseudo inverse Phi which is not the identity
    complexity <- (apply(F,1,function(x) sum(x^2)))^2/apply(F,1,function(x)sum(x^4))   #added  05/25/24
    colnames(sm) <- paste0("F",1:nfactors,"*")
    if(!is.null(Phi)) { result <- list(sl = cbind(gprimaryload, sm,h2, u2,p =g.percent), orthog = orth.load, oblique=fload,
        phi =factr, gloading = gload,S.Phi = Phi.S,Call=cl)} else{
    result <- list(sl = cbind(gprimaryload, sm,h2, u2,p=g.percent, com=complexity), orthog = orth.load, oblique=fload,
        phi =factr, gloading = gload,dof=fact$dof,objective=fact$criteria[1],STATISTIC=fact$STATISTIC,PVAL=fact$PVAL,RMSEA=fact$RMSEA,BIC=fact$BIC,rms = fact$rms,crms=fact$crms,n.obs=n.obs,scores=fact$scores,S.Phi = Phi.S,Call=cl )}
    class(result) <- c("psych" ,"schmid")
    return(result)
}


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
     				               if(inherits(ob, as.character("try-error")))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
     				               ob <- Promax(loadings)}
     				                 
     				loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}
  		                             }
 return(list(loadings=loadings,Phi=Phi))
 }
 
 
 #direct Schmid Leiman adapted from Waller (2017)
 directSl <- function(m,nfactors=3,fm="minres",rotate="oblimin",cut=.3){
 cl <- match.call()
 nvar <- ncol(m)
 if(isCorrelation(m)) {C <- m} else { C <- cor(m,use="pairwise")}
 f <- fa(C,nfactors=nfactors,fm=fm,rotate ='none')   #unrotated solution
 #construct the target from the rotated solution
 f.obl <- oblique.rotations(rotate=rotate,loadings = f$loadings)$loadings
 targ <- factor2cluster(f.obl,cut=cut)
 #Waller adjustments to target and factor model
 targ <- cbind(g=rep(1,nvar),targ)
 f0 <- cbind(rep(0,nvar),f$loadings)
 direct <- Procrustes(f0,targ)$loadings   #The Waller Procrustes solution

 colnames(direct) <- c("g",paste0("F",1:nfactors,"*"))  #put some labels in
 class(direct) <- "loadings"
 results <- list(direct=direct,C=C,f=f,targ=targ,Call=cl)
 class(results) <- c("psych","direct")
 return(results)
 }
 
 omegaDirect <- function(m,nfactors=3,fm="minres",rotate="oblimin",cut=.3,plot=TRUE,main="Direct Schmid Leiman"){
  cl <- match.call()
 dsl  <- directSl(m=m,nfactors=nfactors,fm=fm,rotate=rotate,cut=cut)
 direct <- dsl$direct
 m <- dsl$C
 f <- dsl$f
 targ <- dsl$targ
  if(isCorrelation(m)) {C <- m} else { C <- cor(m,use="pairwise")}
 
 sum.g <- sum(direct[,1])
 Vt <- sum(C)  #the total variance in the matrix
 omega.g <-sum.g^2/Vt
 h2 <- rowSums(direct^2) 
 H2 <- sum(h2)
 u2 <-1 - h2
 U2 <- sum(u2)
 
 om.tot <-1 - U2/Vt
 #find subset omegas
 
  omg <- omgo <- omt<-  rep(NA,nfactors+1)
     sub <- apply(direct,1,function(x) which.max(abs(x[2:(nfactors+1)])))
     grs <- 0
     for(group in( 1:nfactors)) {
     groupi <- which(sub==group)
     if(length(groupi) > 0) {
      Vgr <- sum(C[groupi,groupi])
      gr <- sum(direct[groupi,(group+1)])
      grs <- grs + gr^2
      omg[group+1] <- gr^2/Vgr
      omgo[group+1] <- sum(direct[groupi,1])^2/Vgr
      omt[group+1] <- (gr^2+ sum(direct[groupi,1])^2)/Vgr
     }}
     omgo[1] <- sum(direct[,1])^2/sum(C)  #omega h
     omg[1] <- grs/sum(C)  #omega of subscales
     omt[1] <- om.tot 
     om.group <- data.frame(total=omt,general=omgo,group=omg)
     rownames(om.group) <- colnames(direct)[1:(nfactors+1)]
 
 result <- list(loadings=direct,omega.g=omega.g,om.group=om.group,orth.f = f,Target=targ,Call=cl)
 class(result) <-  c("psych" ,"omegaDirect")
if(plot)  omega.diagram(result,sort=TRUE,simple=FALSE,cut=cut,main=main)
 return(result)

 }


