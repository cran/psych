"factor.minres" <- 
function(r,nfactors=1,residuals=FALSE,rotate="varimax",n.obs = NA,scores=FALSE,SMC=TRUE,missing=FALSE,impute="median", min.err = .001,digits=2,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres") {
 cl <- match.call()
  .Deprecated("fa",msg="factor.minres is deprecated.  Please use the fa function.")
 ##first some functions that are internal to factor.minres
  #this does the ULS fitting  
  "fit.residuals.ols" <- function(Psi,S,nf) {
              diag(S) <- 1- Psi
              eigens <- eigen(S)
              eigens$values[eigens$values  < .Machine$double.eps] <- 100 * .Machine$double.eps
       
         if(nf >1 ) {loadings <- eigens$vectors[,1:nf] %*% diag(sqrt(eigens$values[1:nf])) } else {loadings <- eigens$vectors[,1] * sqrt(eigens$values[1] ) }
         model <- loadings %*% t(loadings)
         residual <- (S - model)^2
         diag(residual) <- 0
         error <- sum(residual)
         }
   "fit.residuals.min.res" <- function(Psi,S,nf) {
              diag(S) <-1- Psi
           
              eigens <- eigen(S)
        #loadings <- eigen.loadings(eigens)[,1:nf]
         if(nf >1 ) {loadings <- eigens$vectors[,1:nf] %*% diag(sqrt(eigens$values[1:nf])) } else {loadings <- eigens$vectors[,1] * sqrt(eigens$values[1] ) }
         model <- loadings %*% t(loadings)
         residual <- (S - model)
         diag(residual) <- 0
        
         error <- det(residual)
         }
 #this code is taken (with minor modification to make ULS) from factanal        
 #it does the iterative calls to fit.residuals   
     "min.res" <- function(S,nf) {
          S.smc <- smc(S)
          if(sum(S.smc) == nf) {start <- rep(.5,nf)}  else {start <- 1- S.smc}       
          res <- optim(start, fit.residuals.ols, method = "L-BFGS-B", lower = .005, 
        upper = 1, control = c(list(fnscale = 1, parscale = rep(0.01, 
            length(start)))), nf= nf, S=S )
    Lambda <- FAout(res$par, S, nf)
    result <- list(loadings=Lambda,res=res)
    }
          
 #these  were also taken from factanal        
     FAout <- function(Psi, S, q) {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE)
        L <- E$vectors[, 1L:q, drop = FALSE]
        load <- L %*% diag(sqrt(pmax(E$values[1L:q] - 1, 0)), 
            q)
        diag(sqrt(Psi)) %*% load
    }
    FAfn <- function(Psi, S, q) {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
        e <- E$values[-(1L:q)]
        e <- sum(log(e) - e) - q + nrow(S)
        -e
    }
 
 ## now start the main function
 
 if((fm !="pa") & (fm != "minres")) {message("factor method not specified correctly, minimum residual used  used")
   fm <- "minres" }
 
    n <- dim(r)[2]
    if (n!=dim(r)[1]) {  n.obs <- dim(r)[1]
               if(scores) {x.matrix <- r
             if(missing) {     #impute values 
        miss <- which(is.na(x.matrix),arr.ind=TRUE)
        if(impute=="mean") {
       item.means <- colMeans(x.matrix,na.rm=TRUE)   #replace missing values with means
       x.matrix[miss]<- item.means[miss[,2]]} else {
       item.med   <- apply(x.matrix,2,median,na.rm=TRUE) #replace missing with medians
        x.matrix[miss]<- item.med[miss[,2]]}
        }}
    		r <- cor(r,use="pairwise") # if given a rectangular matrix, then find the correlations first
           } else {
     				if(!is.matrix(r)) {  r <- as.matrix(r)}
     				 sds <- sqrt(diag(r))    #convert covariance matrices to correlation matrices
                     r <- r/(sds %o% sds)  } #added June 9, 2008
    if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0)}

   
   
    r.mat <- r
    Phi <- NULL 
    colnames(r.mat) <- rownames(r.mat) <- colnames(r)
     if(SMC) { 
      if(nfactors < n/2)   {diag(r.mat) <- smc(r) }  else {if (warnings) message("too many factors requested for this number of variables to use SMC, 1s used instead")}  }
    orig <- diag(r)
   
   
    comm <- sum(diag(r.mat))
    err <- comm
     i <- 1
    comm.list <- list()
    if(fm=="pa") {
    while(err > min.err)    #iteratively replace the diagonal with our revised communality estimate
      {
        eigens <- eigen(r.mat,symmetric=symmetric)
        #loadings <- eigen.loadings(eigens)[,1:nfactors]
         if(nfactors >1 ) {loadings <- eigens$vectors[,1:nfactors] %*% diag(sqrt(eigens$values[1:nfactors])) } else {loadings <- eigens$vectors[,1] * sqrt(eigens$values[1] ) }
         model <- loadings %*% t(loadings)
         
         new <- diag(model)
         
         comm1 <- sum(new)
         diag(r.mat) <- new
         err <- abs(comm-comm1)
         if(is.na(err)) {warning("imaginary eigen value condition encountered in fa,\n Try again with SMC=FALSE \n exiting fa")
          break}
         comm <- comm1
         comm.list[[i]] <- comm1
         i <- i + 1
         if(i > max.iter) {if(warnings)  {message("maximum iteration exceeded")}
                err <-0 }
       }
       
       }
       
       if(fm == "minres") { #added April 12, 2009 to do ULS fits
       uls <- min.res(r,nfactors)
       eigens <- eigen(r)  #used for the summary stats
       result$par <- uls$res
       loadings <- uls$loadings
                            }
       
       # a weird condition that happens with the Eysenck data
       #making the matrix symmetric solves this problem
       if(!is.double(loadings)) {warning('the matrix has produced imaginary results -- proceed with caution')
       loadings <- matrix(as.double(loadings),ncol=nfactors) } 
       #make each vector signed so that the maximum loading is positive  - probably should do after rotation
       #Alternatively, flip to make the colSums of loading positive
   if (FALSE) {
   if (nfactors >1) {
    		maxabs <- apply(apply(loadings,2,abs),2,which.max)
   			sign.max <- vector(mode="numeric",length=nfactors)
    		for (i in 1: nfactors) {sign.max[i] <- sign(loadings[maxabs[i],i])}
    		loadings <- loadings %*% diag(sign.max)
   		
    	} else {
    		mini <- min(loadings)
   			maxi <- max(loadings)
    		if (abs(mini) > maxi) {loadings <- -loadings }
    		loadings <- as.matrix(loadings)
    		if(fm == "minres") {colnames(loadings) <- "MR1"} else {colnames(loadings) <- "PA1"}
    	} #sign of largest loading is positive
    }
    
    #added January 5, 2009 to flip based upon colSums of loadings
    if (nfactors >1) {sign.tot <- vector(mode="numeric",length=nfactors)
                 sign.tot <- sign(colSums(loadings))
                 loadings <- loadings %*% diag(sign.tot)
     } else { if (sum(loadings) <0) {loadings <- -as.matrix(loadings)} else {loadings <- as.matrix(loadings)}
             colnames(loadings) <- "MR1" }
     
 
    #end addition
    if(fm == "minres") {colnames(loadings) <- paste("MR",1:nfactors,sep='')	} else {colnames(loadings) <- paste("PA",1:nfactors,sep='')}
    rownames(loadings) <- rownames(r)
    loadings[loadings==0.0] <- 10^-15    #added to stop a problem with varimax if loadings are exactly 0
    
    model <- loadings %*% t(loadings)  
   
    f.loadings <- loadings #used to pass them to factor.stats 
    if(rotate != "none") {if (nfactors > 1) {
   
    
   	if (rotate=="varimax" | rotate=="quartimax") { 
   			rotated <- do.call(rotate,list(loadings))
   			loadings <- rotated$loadings
   			 Phi <- NULL} else { 
     			if ((rotate=="promax")|(rotate=="Promax")) {pro <- Promax(loadings)
     			                loadings <- pro$loadings
     			                Phi <- pro$Phi} else {
     			if (rotate == "cluster") {loadings <- varimax(loadings)$loadings           			
								pro <- target.rot(loadings)
     			              	loadings <- pro$loadings
     			                Phi <- pro$Phi} else {
     			                     
     			if (rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax") {
     				if (!requireNamespace('GPArotation')) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
     				Phi <- NULL} else { ob  <- do.call(rotate,list(loadings) )
     				loadings <- ob$loadings
     				 Phi <- ob$Phi}
     		                             }
     	               }}}
     	  
     }}
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
   
   result <- factor.stats(r,loadings,Phi,n.obs)   #do stats as a subroutine common to several functions

    result$communality <- round(diag(model),digits)
    result$uniquenesses <- round(diag(r-model),digits)
    result$values <- round(eigens$values,digits)
    result$loadings <- loadings

    if(!is.null(Phi)) {result$Phi <- Phi}
    if(fm == "pa") result$communality.iterations <- round(unlist(comm.list),digits)
    
    if(scores) {result$scores <- factor.scores(x.matrix,loadings) }
    result$factors <- nfactors 
    result$fn <- "factor.minres"
    result$fm <- fm
    result$Call <- cl
    class(result) <- c("psych", "fa")
    return(result) }
    
    #modified October 30, 2008 to sort the rotated loadings matrix by the eigen values.
 
 