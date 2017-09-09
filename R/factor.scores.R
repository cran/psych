"factor.scores" <- function(x,f,Phi=NULL,method=c("Thurstone","tenBerge","Anderson","Bartlett","Harman","components"),rho=NULL,impute="none") {
#the normal case is f is the structure matrix and Phi is not specified
#Note that the Grice formulas distinguish between Pattern and Structure matrices
#I need to confirm that I am doing this

    if(length(method) > 1) method <- "tenBerge"   #the default
    if(method=="regression") method <- "Thurstone"
    if(method=="tenberge") method <- "tenBerge"
    if(length(class(f)) > 1) { if(class(f)[2] =="irt.fa" ) f <- f$fa  }
    
     if(!is.matrix(f)) {Phi <- f$Phi
     f <- loadings(f)
      if(ncol(f)==1) {method <- "Thurstone"}
      }
     nf <- dim(f)[2]
      if(is.null(Phi)) Phi <- diag(1,nf,nf)
     if(dim(x)[1] == dim(f)[1]) {r <- as.matrix(x)
         square <- TRUE} else { 
          square <- FALSE
         if(!is.null(rho)) {r <- rho } else {
          r <- cor(x,use="pairwise") #find the correlation matrix from the data
      }}
      
      S <- f %*% Phi   #the Structure matrix 
   switch(method,   
    "Thurstone" = { w <- try(solve(r,S),silent=TRUE )  #these are the factor weights (see Grice eq. 5)
     	if(class(w)=="try-error") {message("In factor.scores, the correlation matrix is singular, an approximation is used")
               r <- cor.smooth(r)}
        
      w <- try(solve(r,S),silent=TRUE)
      if(class(w)=="try-error") {message("I was unable to calculate the factor score weights, factor loadings used instead")
               w <- f}
      colnames(w) <- colnames(f)
      rownames(w) <- rownames(f)
       }, 
      
  "tenBerge" = { #Following Grice equation 8 to estimate scores for oblique solutions (with a correction to the second line where r should r.inv
        L <- f %*% matSqrt(Phi)
        r.5 <- invMatSqrt(r)
       
        r <- cor.smooth(r)
        inv.r <- try(solve(r),silent=TRUE)
        if(class(inv.r)== as.character("try-error"))  {warning("The tenBerge based scoring could not invert the correlation matrix, regression scores found instead")
                                                      ev <- eigen(r)
      ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
        r <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
        diag(r)  <- 1
       w <- solve(r,f)}  else {
        C <- r.5 %*% L %*% invMatSqrt(t(L) %*% inv.r %*% L)    #note that this is the correct formula, per Grice personal communication
        w <- r.5 %*% C %*% matSqrt(Phi)}
        colnames(w) <- colnames(f)
        rownames(w) <- rownames(f)
        },

 
       
       "Harman" = { #Grice equation 10 -- 
     #   m <- t(f)  %*% f  #factor intercorrelations 
     m <- f %*% t(S)  #should be this  (the model matrix)  Revised August 31, 2017
     diag(m) <- 1  #Grice does not say this, but it is necessary to make it work!
       inv.m <- solve(m)
     #  w <- f %*%inv.m  
     w <- inv.m %*% f
       }, 
       
       
        
    "Anderson" =  { #scores for orthogonal factor solution will be orthogonal  Grice Eq 7 and 8
    I <- diag(1,nf,nf)
    h2 <-  diag( f %*% Phi %*% t(f))
    U2 <- 1 - h2
    inv.U2 <- diag(1/U2)
    w <- inv.U2 %*% f %*% invMatSqrt(t(f) %*% inv.U2 %*% r %*% inv.U2 %*% f)
    colnames(w) <- colnames(f)
    rownames(w) <- rownames(f)
    },
    
   "Bartlett" = {    #Grice eq 9  # f should be the pattern, not the structure 
    I <- diag(1,nf,nf)
    h2 <-  diag( f %*% Phi %*% t(f))
    U2 <- 1 - h2
    inv.U2 <- diag(1/U2)
    w <- inv.U2 %*% f %*% (solve(t(f) %*% inv.U2 %*% f))
    colnames(w) <- colnames(f)
    rownames(w) <- rownames(f)
    },
    "none" = {w <- NULL},
    
    "components" = {w <- try(solve(r,f),silent=TRUE )    #basically, just do the regression/Thurstone approach for components
                    w <- f }
    )
    
    
    #now find a few fit statistics
    if(is.null(w)) {results <- list(scores=NULL,weights=NULL)} else {
     R2 <- diag(t(w) %*% S)  #this had been   R2 <- diag(t(w) %*% f)   Corrected Sept 1, 2017
     if(any(R2 > 1) || (prod(!is.nan(R2)) <1) || (prod(R2) < 0) ) {#message("The matrix is probably singular -- Factor score estimate results are likely incorrect")
                      R2[abs(R2) > 1] <- NA
                      R2[R2 <= 0] <- NA
                     }
     #if ((max(R2,na.rm=TRUE) > (1 + .Machine$double.eps)) ) {message("The estimated weights for the factor scores are probably incorrect.  Try a different factor extraction method.")}
      r.scores <- cov2cor(t(w) %*% r %*% w) #what actually is this?
     
    
  if(square) {  #that is, if given the correlation matrix
     class(w) <- NULL
     results <- list(scores=NULL,weights=w)
      results$r.scores <- r.scores 
   	  results$R2 <- R2   #this is the multiple R2 of the scores with the factors
     } else {
         missing <- rowSums(is.na(x))
    if(impute !="none") {
       x <- data.matrix(x)
        miss <- which(is.na(x),arr.ind=TRUE)
        if(impute=="mean") {
       		item.means <- colMeans(x,na.rm=TRUE)   #replace missing values with means
       		x[miss]<- item.means[miss[,2]]} else { 
       		item.med   <- apply(x,2,median,na.rm=TRUE) #replace missing with medians
        	x[miss]<- item.med[miss[,2]]}   #this only works if items is a matrix
     }
      
     if(method !="components") {scores <- scale(x) %*% w } else {  #standardize the data before doing the regression if using factors, 
        scores <- x %*% w}       # for components, the data have already been zero centered and, if appropriate, scaled
     results <- list(scores=scores,weights=w)
     results$r.scores <- r.scores
     results$missing <- missing 
   	  results$R2 <- R2   #this is the multiple R2 of the scores with the factors
     }
     }
   
     return(results) }
     #how to treat missing data?  see score.item
    
     
     
     
"matSqrt" <- function(x) {
   e <- eigen(x)
    e$values[e$values < 0] <- .Machine$double.eps
   sqrt.ev <- sqrt(e$values)   #need to put in a check here for postive semi definite
   result <- e$vectors %*% diag(sqrt.ev) %*% t(e$vectors)
   result}
   
   
"invMatSqrt" <- function(x) {
   e <- eigen(x)
   if(is.complex(e$values)) {warning("complex eigen values detected by invMatSqrt, results are suspect")
                 result <- x
      } else {
      
       e$values[e$values < .Machine$double.eps] <- 100 * .Machine$double.eps
   inv.sqrt.ev <- 1/sqrt(e$values)   #need to put in a check here for postive semi definite
   result <- e$vectors %*% diag(inv.sqrt.ev) %*% t(e$vectors) }
   result}
