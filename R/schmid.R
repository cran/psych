#corrected estimate of communality, May 21, 2007
#removed "x" from factanal call, June 14, 2007
#added ability to do 2 factors by treating them with equal loadings Jan 2008
#added use of simplimax rotation June 2008
#corrected the sign of group factors to match the original factors
#modified November, 2014 to allow for covariances in addition to correlations.
#also cleaned up the code to take advantage of switch
"schmid" <-
function (model, nfactors = 3, fm = "minres",  digits=2,rotate="oblimin",n.obs=NA,option="equal",Phi=NULL,covar=FALSE,...) 
{
 cl <- match.call()
#if Phi is not Null, then we have been given a factor matrix, otherwise
#model is a correlation matrix, or if not, the correlation matrix is found
#nfactors is the number of factors to extract
      if(!requireNamespace('GPArotation')) {stop("I am sorry, you need to have the  GPArotation package installed")}
  if(is.null(Phi)) {  #the normal case
      normal.case <- TRUE
      nvar <-dim(model)[2]
      if(dim(model)[1] != dim(model)[2]) {n.obs <- dim(model)[1]
                                        if(covar) { model <- cov(model,use="pairwise")} else {model <- cor(model,use="pairwise") }
                                         }
                                          
     if (fm =="pc") {
        fact <- principal(model, nfactors,n.obs=n.obs,...)
    } else {if ((fm == "pa") |(fm =="minres") | (fm =="wls")  |(fm =="minres") |(fm =="ml")|(fm =="mle")  |(fm =="gls") |(fm =="minchi")) {fact <- fa(model, nfactors,n.obs=n.obs,rotate="varimax",fm=fm,covar=covar) } else {
     
        stop("The method of factor extraction you specified is not available")
        
    }}
     orth.load <- loadings(fact)
    } else {model <- as.matrix(model)
            Phi <- as.matrix(Phi)
            fact <- model %*% Phi  #find the orthogonal matrix from the oblique pattern and the Phi matrix
            orth.load <- fact
            nfactors <- dim(fact)[2]
            normal.case <-FALSE}
   
    colnames(orth.load)  <- paste("F",1:nfactors,sep="")
    if(nfactors == 1) { message("Omega_h for 1 factor is not meaningful, just omega_t")
                        obminfact <-list(loadings= orth.load)
                       factr <- 1
                       
           } else {  #the normal case is nfactors > 2
switch(rotate,
    simplimax = {obminfact <- GPArotation::simplimax(orth.load)},
    promax  =    {obminfact  <- Promax(orth.load)
     								 rotmat <- obminfact$rotmat
                   						Phi <- obminfact$Phi
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
           							        if(class(obminfact)== as.character("try-error")) {obminfact <- Promax(orth.load)   #special case for examples with exactly 2 orthogonal factors
           							        message("\nThe oblimin solution failed, Promax used instead.\n")                   #perhaps no longer necessary with patch to GPForth and GPFoblq in GPArotation
           							        rotmat <- obminfact$rotmat
                   						    Phi <- obminfact$Phi}}
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
                        factr <- Phi} else {
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
      			 message("\nThree factors are required for identification -- general factor loadings set to be equal. \nProceed with caution. \nThink about redoing the analysis with alternative values of the 'option' setting.\n")} else { if(option=="first") {
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
    
    
    colnames(sm) <- paste("F",1:nfactors,"*",sep="")
    if(!is.null(Phi)) { result <- list(sl = cbind(gprimaryload, sm,h2, u2,p =g.percent), orthog = orth.load, oblique=fload,
        phi =factr, gloading = gload,Call=cl)} else{
    result <- list(sl = cbind(gprimaryload, sm,h2, u2,p=g.percent), orthog = orth.load, oblique=fload,
        phi =factr, gloading = gload,dof=fact$dof,objective=fact$criteria[1],STATISTIC=fact$STATISTIC,PVAL=fact$PVAL,RMSEA=fact$RMSEA,BIC=fact$BIC,rms = fact$rms,crms=fact$crms,n.obs=n.obs,scores=fact$scores,Call=cl )}
    class(result) <- c("psych" ,"schmid")
    return(result)
}
