
#A number of estimates of unidimensionality
#Developed March 9. 2017
#Modified August 3 to consider two more estimates
#cleaned up Sept 18, 2023 to match ms
#reversed keys and x to be consistent with other functions

"unidim" <- function(keys=NULL,x=NULL,cor="cor",use="pairwise",fm="minres", correct=.5,
     check.keys=TRUE,n.obs=NA, nfactors=3) {
   cl <- match.call()
    covar <- FALSE 
   if(is.null(n.obs)) {
   
    n.obs <- NROW(x)  
   use <- "pairwise"}
    if(!is.null(x)) {n.keys <- length(keys)
    } else {x <- keys
          keys <- NULL 
           n.keys <- 1}
    all.x <- x
   results <- list()
   fits <- list()

   
   for(scale in 1:n.keys) { if(!is.null(keys)) {
    select <- selectFromKeys(keys[scale]) 
# select <- keys[[scale]]
 #     flipper <- rep(1,length(select))
#          flipper[grep("-",select)] <- -1
#          if(is.numeric(select)) {select <- abs(select) } else {
#          select <- sub("-","",unlist(select)) }
   if(isCorrelation(all.x)) {x <- all.x[select,select] } else {x <- all.x[,select]}
}  else {flipper <- rep(1,ncol(x))} #this allows us to handle multiple scales 
   
 if(!isCorrelation(x) ) {
     if(is.na(n.obs)) n.obs <- NROW(x)
        switch(cor, 
       cor = { x <- cor(x,use=use)
              covar=FALSE},
       cov = {x <- cov(x,use=use) 
              covar <- TRUE},
      
       spearman = {x <- cor(x,use=use,method="spearman")},
       kendall = {x <- cor(x,use=use,method="kendall")},
       tet = {x <- tetrachoric(x,correct=correct)$rho},
       tetrachoric = {x <- tetrachoric(x,correct=correct)$rho},
        poly = {x <- polychoric(x,correct=correct)$rho},
       polychoric = {x <- polychoric(x,correct=correct)$rho},
       mixed = {x <- mixedCor(x,use=use,correct=correct)$rho}
       
       )}
 #X is now a correlation matrix    or a covariance matrix 
  
  
    f1 <- fa(x,n.obs=n.obs,covar=covar,fm=fm)  # a one factor solution
  #But to get a reasonable estimate of ECV we need more factors
  nvar <- NCOL(x)
  if(nvar < 4 ) nfactors<-1
  fn <- fa(x,nfactors,n.obs=n.obs,covar=covar,fm=fm,rotate="none") 
  if(f1$dof < 1 ) f1$RMSEA <- NA
  if(is.null(f1$BIC)) f1$BIC <- NA
  if(is.null(f1$CFI)) f1$CFI <- NA
  if(is.null(f1$TLI)) f1$TLI <- NA
  if(is.null(f1$RMSEA)) f1$RMSEA <- NA
  f1f2 <- fn$e.values[1]/fn$e.values[2]
  fa.stats <- list(dof=f1$dof, chi=f1$chi,RMSEA = f1$RMSEA[1],Vaccounted=f1$Vaccounted[2,1],
         TLI = f1$TLI, BIC = f1$BIC, R2=f1$R2, CFI = f1$CFI, ECV=fn$ECV[1], f1f2) 

  g <- sum(f1$model)  # sum(f1$loadings %*% t(f1$loadings))
  n <- NCOL(x)
  Vt <- sum(x)
  om.g <- g/Vt                          #model/ r
  om.t <- (Vt - sum(f1$uniqueness))/Vt   #total reliability 

 uni.orig <- g/ (Vt - sum(f1$uniqueness))  #raw unidimensionality
 
   
 #now, to find traditional alpha, we need to flip negative items
 if(check.keys | n.keys == 1) { flipper <- rep(1,n)
 flipper[sign(f1$loadings ) < 0] <- -1 }
  x <- diag(flipper) %*% x %*% diag(flipper)
  Vt <- sum(x)
 median.r <- median(x[lower.tri(x)])   #lower tri does not count diagonal
 alpha.std <-  (1- tr(x)/Vt)*(n/(n-1))
  av.r <- (Vt-tr(x))/(n*(n-1))
  omega.flip <- sum(diag(flipper) %*% f1$model %*% diag(flipper))/Vt
  omega.total.flip <-  (Vt - sum(f1$uniqueness))/Vt
  flipperped.loadings <- flipper * f1$loadings
  g.flipperped <- sum(flipperped.loadings%*% t(flipperped.loadings))
  uni.flipper <- g.flipperped/(Vt - sum(f1$uniqueness))
  
# How well does the alpha model predict the correlation matrix?
alpha.res <- sum(lower.tri(x)* (x-av.r)^2)/sum(lower.tri(x) * x^2)   #squared residual /squared original
# uni <- uni.flipper * (1-alpha.res)
uni <- f1$fit.off * (1-alpha.res)
 

  stats <- list(u=uni,alpha.res=1-alpha.res, fit.off= f1$fit.off, 
           alpha=alpha.std,av.r = av.r,median.r = median.r,cfi=f1$CFI,ECV=fn$ECV[1], f1f2,
           uni.flipper = uni.flipper, uni=uni.orig,om.g=om.g, omega.pos = omega.flip,om.t=om.t,
         om.total.flip= omega.total.flip)
  if(!is.null(keys)) {results[[names(keys[scale])]]<- stats 
                      fits[[names(keys[scale])]] <- fa.stats} else {results <- stats
                      fits <- fa.stats}
  }
 if(scale==1) {ncol1 <- length(results) } else {ncol1 <- length(results[[1]])}
  
  temp <- matrix(unlist(results),ncol=ncol1,byrow=TRUE)
 
  colnames(temp) <- c("u","tau","con","alpha","av.r","median.r","CFI","ECV","F1/F2","Unidim.A","Unidim","model","model.A", "total", "total.A")
  rownames(temp) <- names(keys)
 if(scale==1) {ncol1 <- length(fits) } else {ncol1 <- length(fits[[1]])}
  fits <- matrix(unlist(fits),ncol=ncol1,byrow=TRUE)
  colnames(fits) <- c("dof" , "chisq" , "RMSEA" ,"Vaccounted","TLI","BIC", "R2","CFI","ECV","F1/F2")
  rownames(fits) <- names(keys)
  if(check.keys) {
  results <- list(uni=temp[,1:9],fa.stats=fits)} else {results <- list(uni=temp,fa.stats=fits)}
  results$Call <- cl
  class(results) <- c("psych","unidim")
  
  return(results)
  }
  
  print_psych.unidim <- function(x,digits=2) {
  cat("\nA measure of unidimensionality \n Call: ")
  print(x$Call)
  
  cat("\nUnidimensionality index = \n" )
  print(round(x$uni,digits=digits))
  
 cat("\nunidim adjusted index reverses negatively scored items.")
 cat("\nalpha ","  Based upon reverse scoring some items.")
 cat ("\naverage and median  correlations are based upon reversed scored items") 
     }
  
  
  
  