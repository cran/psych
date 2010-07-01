#a function to do principal axis, minres,  weighted least squares and maximimum likelihood factor analysis
#basically, just combining the three separate functions
#the code for wls and minres is adapted from the factanal function 
#the optimization function in ml is taken almost directly from the factanal function
#created May 28, 2009 
#modified June 7, 2009 to add gls fitting
#modified June 24, 2009 to add ml fitting
#modified March 4, 2010 to allow for factoring of covariance matrices rather than correlation matrices
#this itself is straight foward, but the summary stats need to be worked on
"fa" <- 
function(r,nfactors=1,n.obs = NA,rotate="oblimin",scores=FALSE,residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1,...) {
 cl <- match.call()
 control <- NULL   #if you want all the options of mle, then use factanal
 
 ##first some functions that are internal to fa
 #this does the WLS or ULS fitting  depending upon fm 
  "fit.residuals" <- function(Psi,S,nf,S.inv,fm) {
              diag(S) <- 1- Psi
              if(!is.null(S.inv)) sd.inv <- diag(1/diag(S.inv))
              eigens <- eigen(S)
              eigens$values[eigens$values  < .Machine$double.eps] <- 100 * .Machine$double.eps
       
         if(nf >1 ) {loadings <- eigens$vectors[,1:nf] %*% diag(sqrt(eigens$values[1:nf])) } else {loadings <- eigens$vectors[,1] * sqrt(eigens$values[1] ) }
         model <- loadings %*% t(loadings)
       #weighted least squares weights by the importance of each variable   
       if(fm == "wls" ) {residual <- sd.inv %*% (S- model)^2 %*% sd.inv} else {if (fm=="gls") {residual <- (S.inv %*%(S - model))^2 } else {residual <- (S - model)^2 #this last is the uls case
       if(fm == "minres") diag(residual) <- 0   #this is minimum residual factor analysis, ignore the diagonal
       }}  # the uls solution usually seems better?
        # 
         error <- sum(residual)
         }
  
 #this next section is taken (with minor modification to make ULS, WLS or GLS) from factanal        
 #it does the iterative calls to fit.residuals 
 #modified June 7, 2009 to add gls fits
 #Modified December 11, 2009 to use first derivatives from formula rather than emprical.  this seriously improves the speed.
     "fit" <- function(S,nf,fm,covar) {
          S.smc <- smc(S,covar)
           if((fm=="wls") | (fm =="gls") ) {S.inv <- solve(S)} else {S.inv <- NULL}
           if(!covar &&(sum(S.smc) == nf)) {start <- rep(.5,nf)}  else {start <- diag(S)- S.smc}
                    #initial communality estimates are variance - smc 
                    
           if(fm=="ml")  {res <- optim(start, FAfn, FAgr, method = "L-BFGS-B",
                          	lower = .005, upper = 1,
                          	control = c(list(fnscale=1,
                 			parscale = rep(0.01, length(start))), control),
                 			nf = nf, S = S)
                 } else {
                		 	res <- optim(start, fit.residuals,gr=FAgr.minres, method = "L-BFGS-B", lower = .005, 
                  			upper = 1, control = c(list(fnscale = 1, parscale = rep(0.01, 
                  			length(start)))), nf= nf, S=S, S.inv=S.inv,fm=fm )
                  		}
   
   if((fm=="wls") | (fm=="gls") ) {Lambda <- FAout.wls(res$par, S, nf)} else { Lambda <- FAout(res$par, S, nf)}
    result <- list(loadings=Lambda,res=res,S=S)
    }
    
 ## the next two functions are taken directly from the factanal function in order to include maximum likelihood as one of the estimate procedures
 
   FAfn <- function(Psi, S, nf)
    {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
        e <- E$values[-(1:nf)]
        e <- sum(log(e) - e) - nf + nrow(S)
       -e
    }
    FAgr <- function(Psi, S, nf)  #the first derivatives
    {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE)
        L <- E$vectors[, 1:nf, drop = FALSE]
        load <- L %*% diag(sqrt(pmax(E$values[1:nf] - 1, 0)), nf)
        load <- diag(sqrt(Psi)) %*% load
        g <- load %*% t(load) + diag(Psi) - S     # g <- model - data
        diag(g)/Psi^2                             #normalized 
    }
    
     FAgr.minres <- function(Psi, S, nf,S.inv,fm)  #the first derivatives 
    {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE)
        L <- E$vectors[, 1:nf, drop = FALSE]
        load <- L %*% diag(sqrt(pmax(E$values[1:nf] - 1, 0)), nf)
        load <- diag(sqrt(Psi)) %*% load
        g <- load %*% t(load) + diag(Psi) - S     # g <- model - data
        diag(g)/Psi^2                             #normalized 
    }
          
 #this was also taken from factanal        
    FAout <- function(Psi, S, q) {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE)
        L <- E$vectors[, 1L:q, drop = FALSE]
        load <- L %*% diag(sqrt(pmax(E$values[1L:q] - 1, 0)), 
            q)
        diag(sqrt(Psi)) %*% load
    }
#This is modified from factanal -- the difference in the loadings is that these produce orthogonal loadings, but slightly worse fit  
   FAout.wls <-  function(Psi, S, q) {
        diag(S) <- 1- Psi
        E <- eigen(S,symmetric = TRUE)
        L <- E$vectors[,1L:q,drop=FALSE] %*%  diag(sqrt(E$values[1L:q,drop=FALSE]),q)
        return(L)
    } ## now start the main function
 
 if((fm !="pa") & (fm != "wls")  & (fm != "gls") & (fm != "minres") & (fm != "uls")& (fm != "ml")) {message("factor method not specified correctly, minimum residual (unweighted least squares  used")
   fm <- "minres" }
 
    n <- dim(r)[2]
    if (n!=dim(r)[1]) {  n.obs <- dim(r)[1]
               if(scores) {x.matrix <- r
             if(missing) { #impute values 
        x.matrix <- as.matrix(x.matrix)  #the trick for replacing missing works only on matrices
        miss <- which(is.na(x.matrix),arr.ind=TRUE)
        if(impute=="mean") {
       item.means <- colMeans(x.matrix,na.rm=TRUE)   #replace missing values with means
       x.matrix[miss]<- item.means[miss[,2]]} else {
       item.med   <- apply(x.matrix,2,median,na.rm=TRUE) #replace missing with medians
        x.matrix[miss]<- item.med[miss[,2]]}
        }}
    		if(!covar) {r <- cor(r,use="pairwise")} else {r <- cov(r,use="pairwise")} # if given a rectangular matrix, then find the correlation or covariance first
           } else {
     				if(!is.matrix(r)) {  r <- as.matrix(r)}
     				if(!covar) {
     				sds <- sqrt(diag(r))    #convert covariance matrices to correlation matrices
                     r <- r/(sds %o% sds) #if we remove this, then we need to fix the communality estimates
                    }
                    } #added June 9, 2008
    if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)}

   
   
    r.mat <- r
    Phi <- NULL 
    colnames(r.mat) <- rownames(r.mat) <- colnames(r)
     if(SMC) { 
      if(nfactors < n/2)   {diag(r.mat) <- smc(r,covar=covar) }  else {if (warnings) message("In fa, too many factors requested for this number of variables to use SMC for communality estimates, 1s are used instead")}  }
    orig <- diag(r)
   
   
    comm <- sum(diag(r.mat))
    err <- comm
     i <- 1
    comm.list <- list()
    
    #principal axis is an iterative eigen value fitting
    if(fm=="pa") {
   	 	e.values <- eigen(r,symmetric=symmetric)$values   #store the original solution
    	while(err > min.err)    #iteratively replace the diagonal with our revised communality estimate
     	 {
       	 eigens <- eigen(r.mat,symmetric=symmetric)
       	  if(nfactors >1 ) {loadings <- eigens$vectors[,1:nfactors] %*% diag(sqrt(eigens$values[1:nfactors])) } else {loadings <- eigens$vectors[,1] * sqrt(eigens$values[1] ) }
        	 model <- loadings %*% t(loadings)
        	 new <- diag(model)       
         	comm1 <- sum(new)
         	diag(r.mat) <- new
        	 err <- abs(comm-comm1)
        	 if(is.na(err)) {warning("imaginary eigen value condition encountered in factor.pa,\n Try again with SMC=FALSE \n exiting factor.pa")
          break}
        	 comm <- comm1
        	 comm.list[[i]] <- comm1
         	i <- i + 1
         	if(i > max.iter) { 
         	         if(warnings)  {message("maximum iteration exceeded")}
                     err <-0 }
          }  #end of while loop  
          eigens <- eigens$values
       } 
       
       if((fm == "wls") | (fm=="minres") | (fm=="gls") | (fm=="uls")|(fm== "ml")) { 
       uls <- fit(r,nfactors,fm,covar=covar)
       
       e.values <- eigen(r)$values  #eigen values of pc: used for the summary stats --  
       result$par <- uls$res
      
       loadings <- uls$loadings
       model <- loadings %*% t(loadings)
       S <- r
       diag(S) <- diag(model)   #communalities from the factor model 
       eigens <- eigen(S)$values
       
                            }
       
       # a weird condition that happens with poor data
       #making the matrix symmetric solves this problem
       if(!is.real(loadings)) {warning('the matrix has produced imaginary results -- proceed with caution')
       loadings <- matrix(as.real(loadings),ncol=nfactors) } 
       #make each vector signed so that the maximum loading is positive  -  should do after rotation
       #Alternatively, flip to make the colSums of loading positive
   
   
    if (nfactors >1) {sign.tot <- vector(mode="numeric",length=nfactors)
                 sign.tot <- sign(colSums(loadings))
                 sign.tot[sign.tot==0] <- 1
                 loadings <- loadings %*% diag(sign.tot)
     } else { if (sum(loadings) <0) {loadings <- -as.matrix(loadings)} else {loadings <- as.matrix(loadings)}
             colnames(loadings) <- "MR1" }
     
 
 
    if(fm == "wls") {colnames(loadings) <- paste("WLS",1:nfactors,sep='')	} else {if (fm=="pa")  {colnames(loadings) <- paste("PA",1:nfactors,sep='')} else  {if (fm=="gls")  {colnames(loadings) <- paste("GLS",1:nfactors,sep='')} else {if (fm=="ml")  {colnames(loadings) <- paste("ML",1:nfactors,sep='')} else {colnames(loadings) <- paste("MR",1:nfactors,sep='')} }}}
    rownames(loadings) <- rownames(r)
    loadings[loadings==0.0] <- 10^-15    #added to stop a problem with varimax if loadings are exactly 0
   
    model <- loadings %*% t(loadings)  
    
    f.loadings <- loadings #used to pass them to factor.stats 
    if(rotate != "none") {if (nfactors > 1) {
   
    
   	if (rotate=="varimax" |rotate=="Varimax" | rotate=="quartimax" | rotate =="bentlerT" | rotate =="geominT") { 
   	       #varimax is from the stats package, Varimax is from GPArotations
   			rotated <- do.call(rotate,list(loadings,...))
   			loadings <- rotated$loadings
   			 Phi <- NULL} else { 
     			if ((rotate=="promax")|(rotate=="Promax") ) {pro <- Promax(loadings)
     			                loadings <- pro$loadings
     			                Phi <- pro$Phi} else {
     			if (rotate == "cluster") {loadings <- varimax(loadings)$loadings           			
								pro <- target.rot(loadings,...)
     			              	loadings <- pro$loadings
     			                Phi <- pro$Phi} else {
     			                     
     			if (rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax" | rotate =="geominQ"  | rotate =="bentlerQ") {
     				if (!require(GPArotation)) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
     				Phi <- NULL} else { ob  <- try(do.call(rotate,list(loadings,...) ))
     				          if(class(ob)== as.character("try-error"))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
     				                  ob <- Promax(loadings)}
     				                 
     				loadings <- ob$loadings
     				 Phi <- ob$Phi}
     		                             }
     	               }}}
     	  
     }}
    signed <- sign(colSums(loadings))
    signed[signed==0] <- 1
    loadings <- loadings %*% diag(signed)  #flips factors to be in positive direction but loses the colnames
    if(!is.null(Phi)) {Phi <- diag(signed) %*% Phi %*% diag(signed) }  #added October 20, 2009 to correct bug found by Erich Studerus
     if(fm == "wls") {colnames(loadings) <- paste("WLS",1:nfactors,sep='')	} else {if (fm=="pa")  {colnames(loadings) <- paste("PA",1:nfactors,sep='')} else  {if (fm=="gls")  {colnames(loadings) <- paste("GLS",1:nfactors,sep='')} else {if (fm=="ml")  {colnames(loadings) <- paste("ML",1:nfactors,sep='')} else {colnames(loadings) <- paste("MR",1:nfactors,sep='')} }}}
        #just in case the rotation changes the order of the factors, sort them
        #added October 30, 2008
       
   if(nfactors >1) {
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated,decreasing=TRUE)
    loadings <- loadings[,ev.order]}
    rownames(loadings) <- colnames(r)
    if(!is.null(Phi)) {Phi <- Phi[ev.order,ev.order] } #January 20, 2009 but, then, we also need to change the order of the rotation matrix!
    class(loadings) <- "loadings"
    if(nfactors < 1) nfactors <- n
   
   result <- factor.stats(r,loadings,Phi,n.obs=n.obs,alpha=alpha)   #do stats as a subroutine common to several functions

    result$communality <- diag(model)
    result$uniquenesses <- diag(r-model)
    result$values <-  eigens
    result$e.values <- e.values  
    result$loadings <- loadings
    result$fm <- fm  #remember what kind of analysis we did
    
    if(!is.null(Phi)) {result$Phi <- Phi}
    if(fm == "pa") result$communality.iterations <- unlist(comm.list)
    
    if(scores) {result$scores <- factor.scores(x.matrix,loadings) }
    result$factors <- nfactors 
    result$fn <- "fa"
    result$fm <- fm
    result$Call <- cl
    class(result) <- c("psych", "fa")
    return(result) }
    
    #modified October 30, 2008 to sort the rotated loadings matrix by the eigen values.
    #modified Spring, 2009 to add multiple ways of doing factor analysis
    #corrected, August, 2009 to count the diagonal when doing GLS or WLS - this mainly affects (improves) the chi square
 