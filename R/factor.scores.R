"factor.scores" <- function(x,f,Phi=NULL,method=c("Thurstone","tenBerge","Anderson","Bartlett","Harman","components"),rho=NULL) {
    if(length(method) > 1) method <- "tenBerge"   #the default
    if(method=="regression") method <- "Thurstone"
     if(!is.matrix(f)) {Phi <- f$Phi
     f <- loadings(f)
      }
     nf <- dim(f)[2]
      if(is.null(Phi)) Phi <- diag(1,nf,nf)
     if(dim(x)[1] == dim(f)[1]) {r <- as.matrix(x)
         square <- TRUE} else { 
          square <- FALSE
         if(!is.null(rho)) {r <- rho } else {
          r <- cor(x,use="pairwise") #find the correlation matrix from the data
      }}
       
   switch(method,   
    "Thurstone" = { w <- try(solve(r,f),silent=TRUE )  #these are the factor weights
     if(class(w)=="try-error") {message("In factor.scores, the correlation matrix is singular, an approximation is used")
               r <- cor.smooth(r)}
        
      w <- try(solve(r,f),silent=TRUE)
      if(class(w)=="try-error") {message("I was unable to calculate the factor score weights, factor loadings used instead")
               w <- f}
      colnames(w) <- colnames(f)
       rownames(w) <- rownames(f)

       }, 
      
    "tenBerge" = { #Following Grice equation 8 to estimate scores for oblique solutions
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
       C <- r.5 %*% L %*% invMatSqrt(t(L) %*% inv.r %*% L)
       w <- r.5 %*% C %*% matSqrt(Phi)}
       colnames(w) <- colnames(f)
       rownames(w) <- rownames(f)
       },
       
    "Harman" = { #Grice equation 10 -- 
       m <- t(f)  %*% f  #factor intercorrelations
       inv.m <- solve(m)
       w <- f %*%inv.m  
       },
       
    "Anderson" =  { #scores for orthogonal factor solution will be orthogonal
    I <- diag(1,nf,nf)
    h2 <-  diag( f %*% Phi %*% t(f))
    U2 <- 1 - h2
    inv.U2 <- diag(1/U2)
    w <- inv.U2 %*% f %*% invMatSqrt(t(f) %*% inv.U2 %*% r %*% inv.U2 %*% f)
    colnames(w) <- colnames(f)
    rownames(w) <- rownames(f)
    },
    
   "Bartlett" = {
    I <- diag(1,nf,nf)
    h2 <-  diag( f %*% Phi %*% t(f))
    U2 <- 1 - h2
    inv.U2 <- diag(1/U2)
    w <- inv.U2 %*% f %*% (solve(t(f) %*% inv.U2 %*% f))
    colnames(w) <- colnames(f)
    rownames(w) <- rownames(f)
    },
    "none" = {w <- NULL},
    
    "components" = {
    w <- f }
    )
    
    
    
    #now find a few fit statistics
    if(is.null(w)) {results <- list(scores=NULL,weights=NULL)} else {
     R2 <- diag(t(w) %*% f)
     if(any(R2>1) || (prod(!is.nan(R2)) <1) || (prod(R2) < 0) ) {#message("The matrix is probably singular -- Factor score estimate results are likely incorrect")
                      R2[abs(R2) > 1] <- NA
                      R2[R2 <= 0] <- NA
                     }
     #if ((max(R2,na.rm=TRUE) > (1 + .Machine$double.eps)) ) {message("The estimated weights for the factor scores are probably incorrect.  Try a different factor extraction method.")}
      r.scores <- cov2cor(t(w) %*% r %*% w) 
     
    
  if(square) {
     class(w) <- NULL
     results <- list(scores=NULL,weights=w)
      results$r.scores <- r.scores 
   	  results$R2 <- R2   #this is the multiple R2 of the scores with the factors
     } else {
      
     scores <- scale(x) %*% w    #standardize the data before doing the regression
     results <- list(scores=scores,weights=w)
     results$r.scores <- r.scores 
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
