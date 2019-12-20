#created Sept 4, 2017 to try to do exploratory sem by factoring sets 1 and 2
#and then linking the two sets
#slightly improved December 15, 2018 to label the factors better

"esem" <- function(r,varsX,varsY,nfX=1,nfY=1,n.obs=NULL,fm="minres",rotate="oblimin",plot=TRUE,cor="cor",use="pairwise",weight=NULL,...) {
if(is.null(colnames(r))) colnames(r) <- rownames(r) <- paste0("V",1:ncol(r))

#vars <- order(c(vars1,vars2))
 cl <- match.call()
 if(is.numeric(varsX)) varsX <- colnames(r)[varsX]
 if(is.numeric(varsY)) varsY <- colnames(r)[varsY] 
vars <- c(varsX,varsY)
if(is.null(n.obs)) n.obs <- NA

if(!isCorrelation(r)) {#find the correlations
n.obs <- nrow(r) 
r <- r[,vars]  #we organize the data to be just the ones we want, that is, ignore variables not included in the model
   switch(cor, 
       cor = {r <- cor(r,use=use)},
       cov = {r <- cov(r,use=use) 
              covar <- TRUE},
       wtd = { r <- cor.wt(r,w=weight)$r},
       tet = {r <- tetrachoric(r)$rho},
       poly = {r <- polychoric(r)$rho},
       mixed = {r <- mixed.cor(r,use=use)$rho},
       Yuleb = {r <- YuleCor(r,,bonett=TRUE)$rho},
       YuleQ = {r <- YuleCor(r,1)$rho},
       YuleY = {r <- YuleCor(r,.5)$rho } 
       )
       }
#varnames <- colnames(r)[vars]
varnames <- vars
rownames(r) <- colnames(r)
R <- r[varnames,varnames]  #This reorganizes R so that it is the order of the selected variables
 nX <- length(varsX)
 nY <-length(varsY)
 df1 <- nX *(nX-1)/2 - nfX * nX + nfX * (nfX-1)/2
  df2 <- nY *( nY-1)/2 - nfY * nY + nfY * (nfY-1)/2 

  f1 <- fa.extend(R,nfX,ov=varsX,ev=varsY,fm=fm,rotate=rotate,...)
 
  loads1 <- f1$loadings[varnames,,drop=FALSE]
  
  
  S1 <- f1$Structure[varnames,,drop=FALSE]
 if(!is.null(ncol(S1)))  colnames(loads1) <- colnames(S1) <- paste0("X",1:ncol(loads1))
  Phi1 <- f1$Phi
  f2 <- fa.extend(R,nfY,ov=varsY,ev=varsX,fm=fm,rotate=rotate,...)

  loads2  <- f2$loadings[varnames,,drop=FALSE]
  S2 <- f2$Structure[varnames,,drop=FALSE]
   if(!is.null(ncol(S2))) {colnames(loads2) <- colnames(S2) <- paste0("Y",1:ncol(loads2))}
    Phi2 <- f2$Phi
  f12 <- cbind(loads1,loads2)
  S12 <- cbind(S1,S2)
  S12 <- as.matrix(S12)
 Phi <- t(S12) %*% solve(R) %*% S12
 loadsX <- f1$loadings[varsX,,drop=FALSE]
 loadsY <- f2$loadings[varsY,,drop=FALSE]
 colnames(loadsX) <- paste0("X",1:ncol(loadsX)) 
 colnames(loadsY) <- paste0("Y",1:ncol(loadsY)) 

# loadsX <- f1$loadings[colnames(R)[varsX],,drop=FALSE]
# loadsY <- f2$loadings[colnames(R)[varsY],,drop=FALSE]
 diag(Phi) <- 1
 #now, a kludge to make it better -- but not actually, so dropped 
if(FALSE){if(!is.null(Phi1)) Phi[1:nfX,1:nfX] <- Phi1
if(!is.null(Phi2))  Phi[(nfX+1):(nfX+nfY),(nfX+1):(nfX+nfY)] <- Phi2
}
 
  result <- esem.stats(R,f12,S12,Phi,n.obs=n.obs) 
  result$n.obs <- n.obs
  result$loadings <- f12
  result$Structure <- S12 
  result$loadsX <- loadsX
  result$loadsY <- loadsY
  result$PhiX <- Phi1
  result$PhiY <- Phi2
   result$esem.dof <- df1 + df2
  result$fm <- fm
  result$fx <- f1$fo
 result$fy <- f2$fo
 result$Phi <- Phi

 result$Call <- cl
 class(result) <- c("psych","esem")
 if(plot) esem.diagram(result)
 
 return(result)}
 
 
 "esem.stats" <- function(r,f,s,phi,n.obs=NA) {
 r <- as.matrix(r)
 n <- ncol(r)
 nfactors <- ncol(f)
 if(is.null(nfactors)) nfactors <- 1
 Sp <- s %*% solve(phi)  #this is not quite f, but is better
 model <- Sp %*% t(s)  #this works better than  model <- f %*%  t(s)

residual <- r - model
result <- list()
 result$communality <- diag(model)
 result$sumsq <- diag(t(s) %*% (Sp))

 r2 <- sum(r*r)
 rstar2 <- sum(residual*residual)
  #Alternatively, we can recognize that we are not estimating all of these
  #this is results$esem.dof
 result$dof <- dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
 
r2.off <- r2 - tr(r)
diag(residual) <- 0
rstar.off <- sum(residual^2)
 
 #   
#   m.inv.r <- try(solve(model,r),silent=TRUE) #modified Oct 30, 2009 to perhaps increase precision -- #modified 2015/1/2 to use try
#       
#      if(inherits(m.inv.r,"try-error") {warning("the model inverse times the r matrix is singular, replaced with Identity matrix which means fits are wrong")
#            sum}
#     if(is.na(n.obs)) {result$n.obs=NA 
#     			      result$PVAL=NA} else {result$n.obs=n.obs}
 m.inv.r <- diag(1,n,n)   #this is because the m.inv.r is not estimated 
  #  result$dof <-  n * (n-1)/2 - n * nfactors + (nfactors *(nfactors-1)/2)
   
    result$objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n   #this is what Tucker Lewis call F
result$objective <- rstar2   #because the normal way doesn't work
#     if(is.infinite(result$objective)) {result$objective <- rstar2
#                                        message("The determinant of the smoothed correlation was zero.\nThis means the objective function is not defined.\nChi square is based upon observed residuals.")}
    result$criteria <- c("objective"=result$objective,NA,NA)
   
     if (!is.na(n.obs)) {result$STATISTIC <-  chisq <- result$objective * ((n.obs-1) -(2 * n + 5)/6 -(2*nfactors)/3) #from Tucker  and from factanal
    # if (!is.na(n.obs)) {result$STATISTIC <-  chisq <- result$objective * ((n.obs-1)) #from Fox and sem
          if(!is.nan(result$STATISTIC)) if (result$STATISTIC <0) {result$STATISTIC <- 0}  
   			if (result$dof > 0) {result$PVAL <- pchisq(result$STATISTIC, result$dof, lower.tail = FALSE)} else {result$PVAL <- NA}
    		}
      
 result$ENull <- r2.off * n.obs  #the empirical null model
 result$null.dof <- n * (n-1)
  result$chi <- rstar.off * n.obs  #this is the empirical chi square
  result$rms <- sqrt(rstar.off/(n*(n-1)))  #this is the empirical rmsea                      
  result$nh <- n.obs
                             if (result$dof > 0) {result$EPVAL <- pchisq(result$chi,        result$dof, lower.tail = FALSE)
result$crms <- sqrt(rstar.off/(2*result$dof) )
 result$EBIC <- result$chi - result$dof * log(n.obs) 
 result$ESABIC <- result$chi - result$dof * log((n.obs+2)/24) } else {result$EPVAL <- NA
 result$crms <- NA
result$EBIC <- NA
 result$ESABIC <- NA}
                                 
     result$fit <-1-rstar2/r2
    result$fit.off <- 1-rstar.off/r2.off
    result$sd <- sd(as.vector(residual)) #this is the non-sample size weighted root mean square residual
    result$factors <- nfactors
    
    result$complexity <- (apply(f,1,function(x)  sum(x^2)))^2/apply(f,1,function(x)sum(x^4))
  result$residual <- residual
    diag(model) <- diag(r) 
    return(result) 
 }
   
  
  "print.psych.esem" <-function(x,digits=2,short=TRUE,cut=NULL, suppress.warnings=TRUE,...)  {

 cat("Exploratory Structural Equation Modeling  Analysis using method = ",x$fm )
   cat("\nCall: ")
   print(x$Call)
   nitems <- nrow(x$loadings)
   nfactors <- ncol(x$loadings)
  
   
   cat("\nFor the 'X' set:\n")
  x$loadsX <- as.matrix(x$loadsX)
   print(round(x$loadsX,digits=digits))
 if(!short) { if(!is.null(ncol(x$PhiX))){   cat("\nWith factor intercorrelations of \n")
 print(round(x$PhiX,digits=digits)) }
   }
  
   
     cat("\nFor the 'Y' set:\n")
     x$loadsY <- as.matrix(x$loadsY)
   print(round(x$loadsY,digits=digits))
   if(!short) {
    if(!is.null(ncol(x$PhiY))) {  cat("\nWith factor intercorrelations of \n")
    print(round(x$PhiY,digits=digits))
    }
    }
   
  if(!short) {
     cat('\nStandardized  pattern coefficients on the X and Y sets using Factor Extension\n')
     L <- cbind(x$loadings,x$communality,1-x$communality)
     
     colnames(L)[(ncol(L)-1):ncol(L)] <- c("h2","u2")
   #  print(round(x$loadings,digits=digits),round(x$communality,digits=digits))
     print(round(L,digits=digits))
     
      varex <- rbind("SS loadings" =   x$sumsq)
        varex <- rbind(varex, "Proportion Var" =   x$sumsq/nitems)
           varex <- rbind(varex, "Cumulative Var"=  cumsum( x$sumsq/nitems))
                              varex <- rbind(varex, "Cum. factor Var"=  cumsum( x$sumsq/sum( x$sumsq)))
    cat("\n") 
    print(round(varex, digits))
     

   }
   
   
 #     if(!is.null(x$complexity)) cat("\nMean item complexity = ",round(mean(x$complexity),1))     
#        objective <- x$criteria[1]
#   
   cat("\nCorrelations between the X and Y sets.\n")
   print(round(x$Phi,digits=digits))


    if(!is.null(x$null.dof)) {cat("\nThe degrees of freedom for the null model are ",x$null.dof, " and the empirical chi square  function was ",round(x$ENull,digits),...)}
    
  
    cat("\nThe degrees of freedom for the model are",x$dof," and the empirical chi square function was ",round(x$chi,digits),"\n"," with prob < ", signif(x$EPVAL,digits),"\n" ,...) 
     
     
    if(!is.null(x$rms)) {cat("\nThe root mean square of the residuals (RMSR) is ", round(x$rms,digits),"\n") }
    if(!is.null(x$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$crms,digits),"\n",...) }
    
  
     if((!is.null(x$chi)) && (!is.na(x$chi))) {cat(" with the empirical chi square ", round(x$chi,digits), " with prob < ", signif(x$EPVAL,digits),"\n" ,...)  }
   	
   	 if(!is.na(x$n.obs)) {cat("The total number of observations was ",x$n.obs, " with fitted Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits),"\n",...)}
  
     
   	if(!is.null(x$TLI)) cat("\nTucker Lewis Index of factoring reliability = ",round(x$TLI,digits+1))
   	if(!is.null(x$RMSEA)) {cat("\nRMSEA index = ",round(x$RMSEA[1],digits+1), " and the", (1- x$RMSEA[4])*100,"% confidence intervals are ",round(x$RMSEA[2:3],digits+1),...)  }
if(!is.null(x$EBIC)) {cat("\nEmpirical BIC = ",round(x$EBIC,digits))}
if(!is.null(x$ESABIC)) {cat("\nESABIC = ",round(x$ESABIC,digits))}


if(!is.null(x$fit)) cat("\nFit based upon off diagonal values =", round(x$fit.off,digits))
 	
 	if(short)  cat("\nTo see the item loadings for the X and Y sets combined, and the associated fa output, print with  short=FALSE.\n")
 	 
# cat("\nTo	 report the factor analysis of the X and Y sets with their associated statistics, run  fa on the X and Y sets separately.")
  
  }
  
  
  "interbattery" <- function(r, varsX, varsY, nfX = 1, nfY = 1, n.obs = NULL,cor = "cor", use = "pairwise",weight=NULL) {
   cl <- match.call()
vars <- c(varsX,varsY)
if(is.null(n.obs)) n.obs <- NA

if(ncol(r)  < nrow(r)) {#find the correlations
n.obs <- nrow(r) 
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
       }
varnames <- colnames(r)[vars]
R <- r[varnames,varnames]  #This reorganizes R so that it is the order of the selected 
  r12 <- r[varsX,varsY]
  H1 <- r12 %*% t(r12)
  E1  <- eigen(H1)
  W1 <- E1$vectors[,1:nfX,drop=FALSE]
  gamma1 <- sqrt(E1$values[1:nfX,drop=FALSE])
  A1 <- W1 %*% diag(sqrt(gamma1),ncol=nfX)
  W2 <- t(r12) %*% W1 %*% diag(1/gamma1,ncol=nfX)
  A2 <- W2 %*% diag(sqrt(gamma1))
  As <- colSums(sign(A1))
  flip <- diag(sign(As),ncol=nfX)
  A1 <- A1 %*% flip
    As <- colSums(sign(A2))
  flip <- diag(sign(As),ncol=nfX)
  A2 <- A2 %*% flip
  
  colnames(A1) <- colnames(A2) <- paste0("IB",1:ncol(A1))
  rownames(A1) <- rownames(r12)
  return(list(A1=A1,A2 = A2,loadings=rbind(A1,A2),Call=cl))
  }

  
 