#a function to do principal axis, minres,  weighted least squares and maximimum likelihood factor analysis
#basically, just combining the three separate functions
#the code for wls and minres is adapted from the factanal function 
#the optimization function in ml is taken almost directly from the factanal function
#created May 28, 2009 
#modified June 7, 2009 to add gls fitting
#modified June 24, 2009 to add ml fitting
#modified March 4, 2010 to allow for factoring of covariance matrices rather than correlation matrices
#this itself is straight foward, but the summary stats need to be worked on
#modified April 4, 2011 to allow for factor scores of oblique or orthogonal solution
#In May, 2011, fa was added as a wrapper to do iterations, and the original fa function was changed to fac.  The functionality of fa has not changed.
#Revised November, 2012 to add the minchi option for factoring.  This minimizes the sample size weighted residual matrix
#Revised 1/2/14 to add mclapply (multicore) feature.  Increase in speed is 50\% for two cores, but only 63\% for 4 cores or 76\% for 8 cores
#dropped the fisherz transform on loadings and phis
#6/12/14  Added the ability to find tetrachorics, polychorics, or mixed cors.
#15/1/15  Fixed the way we handle missing and imputation to actually work.
#19/1/15 modified calls to rotation functions to meet CRAN specs using nameSpace
"fa" <- 
function(r,nfactors=1,n.obs = NA,n.iter=1,rotate="oblimin",scores="regression", residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",...) {
 cl <- match.call()
  if(dim(r)[1] == dim(r)[2] ) {if(is.na(n.obs) && (n.iter >1)) stop("You must specify the number of subjects if giving a correlation matrix and doing confidence intervals")
                               #  if(!require(MASS)) stop("You must have MASS installed to simulate data from a correlation matrix")
                                 }
  
 f <- fac(r=r,nfactors=nfactors,n.obs=n.obs,rotate=rotate,scores=scores,residuals=residuals,SMC = SMC,covar=covar,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor, ...=...) #call fa with the appropriate parameters
 fl <- f$loadings  #this is the original

# if(!require(parallel)) {message("Parallels is required to do confidence intervals")}

 nvar <- dim(fl)[1]
 
 if(n.iter > 1) {
 if(is.na(n.obs) ) {n.obs <- f$n.obs} 
 replicates <- list()
 rep.rots <- list()
 
replicateslist <- parallel::mclapply(1:n.iter,function(x) {
# replicateslist <- lapply(1:n.iter,function(x) {
 if(dim(r)[1] == dim(r)[2]) {#create data sampled from multivariate normal with observed correlation
                                      mu <- rep(0, nvar)
                                      #X <- mvrnorm(n = n.obs, mu, Sigma = r, tol = 1e-06, empirical = FALSE)
                                      #the next 3 lines replaces mvrnorm (taken from mvrnorm, but without the checks)
                                      eX <- eigen(r)
                                      X <- matrix(rnorm(nvar * n.obs),n.obs)
                                      X <-  t(eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(X))
                            } else {X <- r[sample(n.obs,n.obs,replace=TRUE),]}
  fs <- fac(X,nfactors=nfactors,rotate=rotate,scores="none",SMC = SMC,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor,...=...) #call fa with the appropriate parameters
  if(nfactors == 1) {replicates <- list(loadings=fs$loadings)} else  {
                    t.rot <- target.rot(fs$loadings,fl)
                
                   if(!is.null(fs$Phi)) {  phis <- fs$Phi  # should we rotate the simulated factor  correlations?
                   #we should report the target rotated phis, not the untarget rotated phis 
                     replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(t.rot$Phi)])   #corrected 6/10/15
                    #replicates <- list(loadings=t.rot$loadings,phis=phis[lower.tri(phis)])
                    }  else 
                   {replicates <- list(loadings=t.rot$loadings)}                
  }})
  

replicates <- matrix(unlist(replicateslist),nrow=n.iter,byrow=TRUE)

means <- colMeans(replicates,na.rm=TRUE)
sds <- apply(replicates,2,sd,na.rm=TRUE)

if(length(means) > (nvar * nfactors) ) {
   means.rot <- means[(nvar*nfactors +1):length(means)]
   sds.rot <-      sds[(nvar*nfactors +1):length(means)]  
   ci.rot.lower <- means.rot + qnorm(p/2) * sds.rot
  ci.rot.upper <- means.rot + qnorm(1-p/2) * sds.rot  
   ci.rot <- data.frame(lower=ci.rot.lower,upper=ci.rot.upper)    } else  {
        rep.rots <- NULL
         means.rot <- NULL
         sds.rot <- NULL
         z.rot <- NULL
         ci.rot <- NULL }
   
   means <- matrix(means[1:(nvar*nfactors)],ncol=nfactors)
   sds <- matrix(sds[1:(nvar*nfactors)],ncol=nfactors)
   tci <- abs(means)/sds
    ptci <- 1-pnorm(tci)
    if(!is.null(rep.rots)) {
   tcirot <- abs(means.rot)/sds.rot
   ptcirot <- 1- pnorm(tcirot)} else  {tcirot <- NULL
                                      ptcirot <- NULL}
ci.lower <-  means + qnorm(p/2) * sds
ci.upper <- means + qnorm(1-p/2) * sds

ci <- data.frame(lower = ci.lower,upper=ci.upper)
class(means) <- "loadings"

colnames(means) <- colnames(sds) <- colnames(fl)
rownames(means) <- rownames(sds) <- rownames(fl)
f$cis <- list(means = means,sds = sds,ci = ci,p =2*ptci, means.rot=means.rot,sds.rot=sds.rot,ci.rot=ci.rot,p.rot = ptcirot,Call= cl,replicates=replicates,rep.rots=rep.rots)
results <- f 
 results$Call <- cl
class(results) <- c("psych","fa.ci")
} else {results <- f
        results$Call <- cl
       class(results) <- c("psych","fa")
       }
return(results)

 }
 #written May 1 2011
 #modified May 8, 2014 to make cis an object in f to make sorting easier


#the main function 

"fac" <- 
function(r,nfactors=1,n.obs = NA,rotate="oblimin",scores="tenBerge",residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",...) {
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
    #use switch to clean up the code
    switch(fm,
    wls = {residual <- sd.inv %*% (S- model)^2 %*% sd.inv},
    gls = {residual <- (S.inv %*%(S - model))^2 } ,
    uls = {residual <- (S - model)^2},  
    minres = {residual <- (S - model)^2
            diag(residual) <- 0},
    minchi = {residual <- (S - model)^2   #min chi does a minimum residual analysis, but weights the residuals by their pairwise sample size
            residual <- residual * np.obs
            diag(residual) <- 0
            })
        
   #     #weighted least squares weights by the importance of each variable   
#        if(fm == "wls" ) {residual <- sd.inv %*% (S- model)^2 %*% sd.inv} else {if (fm=="gls") {residual <- (S.inv %*%(S - model))^2 } else {residual <- (S - model)^2 #this last is the uls case
#        if(fm == "minres") {diag(residual) <- 0}   #this is minimum residual factor analysis, ignore the diagonal
#        if(fm=="minchi") {residual <- residual * np.obs
#                          diag(residual) <- 0 }   #min chi does a minimum residual analysis, but weights the residuals by their pairwise sample size
#        }}  # the uls solution usually seems better than wls or gls?
#         # 
         error <- sum(residual)
         }
  
 #this next section is taken (with minor modification to make ULS, WLS or GLS) from factanal        
 #it does the iterative calls to fit.residuals 
 #modified June 7, 2009 to add gls fits
 #Modified December 11, 2009 to use first derivatives from formula rather than empirical.  This seriously improves the speed.
     "fit" <- function(S,nf,fm,covar) {
          S.smc <- smc(S,covar)
           if((fm=="wls") | (fm =="gls") ) {S.inv <- solve(S)} else {S.inv <- NULL}
           if(!covar &&(sum(S.smc) == nf) && (nf > 1)) {start <- rep(.5,nf)}  else {start <- diag(S)- S.smc}
                    #initial communality estimates are variance - smc  unless smc = 1 
                    
           if(fm=="ml" || fm=="mle" )  {res <- optim(start, FAfn, FAgr, method = "L-BFGS-B",
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
    
 ## the next two functions are taken directly from the factanal function in order to include maximum likelihood as one of the estimation procedures
 
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
        if(fm=="minchi") {g <- g*np.obs}
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
    #np.obs <- NULL   #only returned with a value in case of fm="minchi" 
 if (fm == "mle" || fm =="MLE" || fm == "ML" ) fm <- "ml"  #to correct any confusion
 if (!any(fm %in%(c("pa","wls","gls","minres","minchi", "uls","ml","mle") ))) {message("factor method not specified correctly, minimum residual (unweighted least squares  used")
   fm <- "minres" }
 
     x.matrix <- r
    n <- dim(r)[2]
    if (n!=dim(r)[1]) {  matrix.input <- FALSE  #return the correlation matrix in this case
                       n.obs <- dim(r)[1]
     
        if(missing) { #impute values 
        x.matrix <- as.matrix(x.matrix)  #the trick for replacing missing works only on matrices
        miss <- which(is.na(x.matrix),arr.ind=TRUE)
        if(impute=="mean") {
       item.means <- colMeans(x.matrix,na.rm=TRUE)   #replace missing values with means
       x.matrix[miss]<- item.means[miss[,2]]} else {
       item.med   <- apply(x.matrix,2,median,na.rm=TRUE) #replace missing with medians
        x.matrix[miss]<- item.med[miss[,2]]}
        }
    		#if(fm=="minchi") 
    		np.obs <- count.pairwise(r)    #used if we want to do sample size weighting
    		if(covar) {cor <- "cov"}  
    # if given a rectangular matrix, then find the correlation or covariance 
    #multiple ways of find correlations or covariances
    switch(cor, 
       cor = {r <- cor(r,use=use)},
       cov = {r <- cov(r,use=use) 
              covar <- TRUE},
       tet = {r <- tetrachoric(r)$rho},
       poly = {r <- polychoric(r)$rho},
       mixed = {r <- mixed.cor(r,use=use)$rho},
       Yuleb = {r <- YuleCor(r,,bonett=TRUE)$rho},
       YuleQ = {r <- YuleCor(r,1)$rho},
       YuleY = {r <- YuleCor(r,.5)$rho } 
       )
       
    		
    
           } else { matrix.input <- TRUE #don't return the correlation matrix
                   if(fm=="minchi") { 
                       if(is.null(np.obs)) {fm <- "minres"
                                message("factor method minchi does not make sense unless we know the sample size, minres used instead")
                            }
                   }
                    if(is.na(n.obs) && !is.null(np.obs))         n.obs <- max(as.vector(np.obs))
     				if(!is.matrix(r)) {  r <- as.matrix(r)}
     				if(!covar) {
     				r <- cov2cor(r)  #probably better to do it this way (11/22/2010)
     				#sds <- sqrt(diag(r))    #convert covariance matrices to correlation matrices
                    # r <- r/(sds %o% sds) #if we remove this, then we need to fix the communality estimates
                    }
                    } #added June 9, 2008
                    #does this next line actually do anything?
    if (!residuals) { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,np.obs=np.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),fit=0)} else { result <- list(values=c(rep(0,n)),rotation=rotate,n.obs=n.obs,np.obs=np.obs,communality=c(rep(0,n)),loadings=matrix(rep(0,n*n),ncol=n),residual=matrix(rep(0,n*n),ncol=n),fit=0,r=r)}
    
   
    if(is.null(SMC)) SMC=TRUE   #if we don't specify it, make it true
    r.mat <- r
    Phi <- NULL 
    colnames(r.mat) <- rownames(r.mat) <- colnames(r)
    if(any(is.na(r))) {
       bad <- TRUE
      tempr <-r
      wcl <-NULL
     while(bad) {
     	wc <- table(which(is.na(tempr), arr.ind=TRUE))  #find the correlations that are NA
    	wcl <- c(wcl,as.numeric(names(which(wc==max(wc)))))
    	tempr <- r[-wcl,-wcl]
    	if(any(is.na(tempr))) {bad <- TRUE} else {bad <- FALSE}
         }

     	cat('\nLikely variables with missing values are ',colnames(r)[wcl],' \n')
      	stop("I am sorry: missing values (NAs) in the correlation matrix do not allow me to continue.\nPlease drop those variables and try again." )
       }
     if(is.logical(SMC) )  {
                  if(SMC) {if(nfactors <= n)   {#changed to <= n instead of < n/2 This warning seems to confuse people unnecessarily
                           diag(r.mat) <- smc(r,covar=covar) 
                           }  else {if (warnings) {
                           message("In fa, too many factors requested for this number of variables to use SMC for communality estimates, 1s are used instead")}
                            }   } else { diag(r.mat) <- 1
                }
              } else { diag(r.mat) <- SMC} 
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
        	 if(is.na(err)) {warning("imaginary eigen value condition encountered in fa\n Try again with SMC=FALSE \n exiting fa")
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
       
       if((fm == "wls") | (fm=="minres") |(fm=="minchi") | (fm=="gls") | (fm=="uls")|(fm== "ml")|(fm== "mle")) { 
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
       if(!is.double(loadings)) {warning('the matrix has produced imaginary results -- proceed with caution')
       loadings <- matrix(as.double(loadings),ncol=nfactors) } 
       #make each vector signed so that the maximum loading is positive  -  should do after rotation
       #Alternatively, flip to make the colSums of loading positive
   
   
    if (nfactors >1) {sign.tot <- vector(mode="numeric",length=nfactors)
                 sign.tot <- sign(colSums(loadings))
                 sign.tot[sign.tot==0] <- 1
                 loadings <- loadings %*% diag(sign.tot)
     } else { if (sum(loadings) <0) {loadings <- -as.matrix(loadings)} else {loadings <- as.matrix(loadings)}
             colnames(loadings) <- "MR1" }
     
    
    switch(fm, 
    wls={colnames(loadings) <- paste("WLS",1:nfactors,sep='')	},
    pa= {colnames(loadings) <- paste("PA",1:nfactors,sep='')} ,
    gls = {colnames(loadings) <- paste("GLS",1:nfactors,sep='')},
    ml = {colnames(loadings) <- paste("ML",1:nfactors,sep='')}, 
    minres = {colnames(loadings) <- paste("MR",1:nfactors,sep='')},
    minchi = {colnames(loadings) <- paste("MC",1:nfactors,sep='')})
    
    rownames(loadings) <- rownames(r)
    loadings[loadings==0.0] <- 10^-15    #added to stop a problem with varimax if loadings are exactly 0
   
    model <- loadings %*% t(loadings)  
    
    f.loadings <- loadings #used to pass them to factor.stats 
    
    rot.mat <- NULL
    if(rotate != "none") {if (nfactors > 1) {

if (rotate=="varimax" |rotate=="Varimax" | rotate=="quartimax" | rotate =="bentlerT" | rotate =="geominT" | rotate =="targetT" | rotate =="bifactor"   | rotate =="TargetT"|
                       rotate =="equamax"| rotate =="varimin"|rotate =="specialT" | rotate =="Promax"  | rotate =="promax"| rotate =="cluster" |rotate == "biquartimin" |rotate == "TargetQ"  |rotate =="specialQ" ) {
Phi <- NULL 
switch(rotate,  #The orthogonal cases  for GPArotation + ones developed for psych
  varimax = {rotated <- stats::varimax(loadings)  #varimax is from stats, the others are from GPArotation 
   			         loadings <- rotated$loadings
   			         rot.mat <- rotated$rotmat},
   Varimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	       #varimax is from the stats package, Varimax is from GPArotations
   			#rotated <- do.call(rotate,list(loadings,...))
   			#rotated <- do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...))
   			rotated <- GPArotation::Varimax(loadings)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	quartimax = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	      
   			#rotated <- do.call(rotate,list(loadings))
   			rotated <- GPArotation::quartimax(loadings)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	bentlerT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	       
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::bentlerT(loadings)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	geominT	= {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	      
   			#rotated <- do.call(rotate,list(loadings,...))
   			rotated <- GPArotation::geominT(loadings)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   	targetT = {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   			rotated <- GPArotation::targetT(loadings,Tmat=diag(ncol(loadings)),...)
   			loadings <- rotated$loadings
   			 rot.mat <- t(solve(rotated$Th))} ,
   			
   	 bifactor = {rot <- bifactor(loadings)
   	             loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},  
   	 TargetT =  {if (!requireNamespace('GPArotation')) {stop("I am sorry, to do this rotation requires the GPArotation package to be installed")}
   	            rot <- GPArotation::targetT(loadings,Tmat=diag(ncol(loadings)),...)
   	              loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},
   	equamax =  {rot <- equamax(loadings)
   	              loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))}, 
   	varimin = {rot <- varimin(loadings)
   	            loadings <- rot$loadings
   	            rot.mat <- t(solve(rot$Th))},
   	specialT =  {rot <- specialT(loadings)
   	              loadings <- rot$loadings
   	              rot.mat <- t(solve(rot$Th))}, 
   	Promax =   {pro <- Promax(loadings)
     			loadings <- pro$loadings
     			 Phi <- pro$Phi 
     			 rot.mat <- pro$rotmat},
     promax =   {pro <- stats::promax(loadings)   #from stats
     			 loadings <- pro$loadings
     			  rot.mat <- pro$rotmat
     			  ui <- solve(rot.mat)
     			  Phi <-  cov2cor(ui %*% t(ui))},	
     cluster = 	 {loadings <- varimax(loadings)$loadings           			
								pro <- target.rot(loadings)
     			              	loadings <- pro$loadings
     			                Phi <- pro$Phi
     			                 rot.mat <- pro$rotmat},
     biquartimin =    {ob <- biquartimin(loadings,)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				 rot.mat <- t(solve(ob$Th))}, 
     TargetQ  =  {ob <- TargetQ(loadings,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}, 
     specialQ = {ob <- specialQ(loadings,...)
                    loadings <- ob$loadings
     				 Phi <- ob$Phi
     				 rot.mat <- t(solve(pro$Th))})
     } else {
     #The following oblique cases all use GPArotation			                
     if (rotate =="oblimin"| rotate=="quartimin" | rotate== "simplimax" | rotate =="geominQ"  | rotate =="bentlerQ"  |rotate == "targetQ"  ) {
     				if (!requireNamespace('GPArotation')) {warning("I am sorry, to do these rotations requires the GPArotation package to be installed")
     				    Phi <- NULL} else { 
     				      
     				             ob <- try(do.call(getFromNamespace(rotate,'GPArotation'),list(loadings,...)))
     				               if(class(ob)== as.character("try-error"))  {warning("The requested transformaton failed, Promax was used instead as an oblique transformation")
     				               ob <- Promax(loadings)}
     				                 
     				loadings <- ob$loadings
     				 Phi <- ob$Phi
     				  rot.mat <- t(solve(ob$Th))}
     		                             } else {message("Specified rotation not found, rotate='none' used")}
     	 }
     	} 
     	 }
     	 		
    signed <- sign(colSums(loadings))
    signed[signed==0] <- 1
    loadings <- loadings %*% diag(signed)  #flips factors to be in positive direction but loses the colnames
    if(!is.null(Phi)) {Phi <- diag(signed) %*% Phi %*% diag(signed) }  #added October 20, 2009 to correct bug found by Erich Studerus
  
    switch(fm, 
    wls={colnames(loadings) <- paste("WLS",1:nfactors,sep='')	},
    pa= {colnames(loadings) <- paste("PA",1:nfactors,sep='')} ,
    gls = {colnames(loadings) <- paste("GLS",1:nfactors,sep='')},
    ml = {colnames(loadings) <- paste("ML",1:nfactors,sep='')}, 
    minres = {colnames(loadings) <- paste("MR",1:nfactors,sep='')},
    uls =  {colnames(loadings) <- paste("ULS",1:nfactors,sep='')},
    minchi = {colnames(loadings) <- paste("MC",1:nfactors,sep='')})
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
    if(max(abs(loadings) > 1.0) && !covar) warning(' A Heywood case was detected.  Examine the loadings carefully.') 
    result <- factor.stats(r,loadings,Phi,n.obs=n.obs,np.obs=np.obs,alpha=alpha)   #do stats as a subroutine common to several functions
    result$rotation <- rotate
    result$communality <- diag(model)
    result$uniquenesses <- diag(r-model)
    result$values <-  eigens
    result$e.values <- e.values  
    result$loadings <- loadings
    result$fm <- fm  #remember what kind of analysis we did
    result$rot.mat <- rot.mat
    if(!is.null(Phi) ) {result$Phi <- Phi      #the if statement was incorrectly including oblique.scores.  Fixed Feb, 2012 following a report by Jessica Jaynes
                       Structure <- loadings %*% Phi} else {Structure <- loadings}
                       class(Structure) <- "loadings"
                       result$Structure <- Structure #added December 12, 2011   
                      
    if(fm == "pa") result$communality.iterations <- unlist(comm.list)
   
    if(oblique.scores) {result$scores <- factor.scores(x.matrix,f=loadings,Phi=Phi,method=scores) } else {result$scores <- factor.scores(x.matrix,f=Structure,method=scores)}

    result$weights <- result$scores$weights
    result$scores <- result$scores$scores
        if(!is.null(result$scores)) colnames(result$scores) <- colnames(loadings) #added Sept 27, 2013
    result$factors <- nfactors 
    result$r <- r   #save the correlation matrix 
    result$np.obs <- np.obs
    result$fn <- "fa"
    result$fm <- fm
    result$Call <- cl
    class(result) <- c("psych", "fa")
    return(result) }
    
    #modified October 30, 2008 to sort the rotated loadings matrix by the eigen values.
    #modified Spring, 2009 to add multiple ways of doing factor analysis
    #corrected, August, 2009 to count the diagonal when doing GLS or WLS - this mainly affects (improves) the chi square
    #modified April 4, 2011 to find the factor scores of the oblique factors
    #modified December 12, 2011 to report structure coefficients as well as pattern (loadings)
   #modified February 11, 2013 to correctly treat SMC=FALSE as 1s instead of 0s.
   #modified spring, 2015 to use switch in the rotation options
   #modified August 25, 2015 to add rot.mat as output