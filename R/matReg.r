
#a helper function to find regressions from covariances
#May 6, 2016
#Fixed November 29, 2018 to handle se of partialed variables correctly
#modified September 25, 2019 to find intercepts and standard errors using momements
#modified even more November 25, 2019 to return the intercepts and R2
matReg <- function(x,y=NULL,C=NULL,m=NULL,z=NULL,n.obs=0,means=NULL,std=FALSE,raw=TRUE,part=FALSE) {
#first, allow for formula based calls
 cl <- match.call()
    #convert names to locations if given formula input
    prod <- ex <-  NULL   #in case we do not have formula input
   #first, see if they are in formula mode  
   if(inherits(x,"formula")) {
   ps <- fparse(x)
   y <- ps$y          #[c(y, x, z, ex)]
   x <- ps$x
   data <- C  #to make the following easier to adapt from setCor
   med <- ps$m #but, mediation is not done here, so we just add this to x
  # if(!is.null(med)) x <- c(x,med)   #not  necessary, because we automatically put this in
   prod <- ps$prod
   z <- ps$z   #do we have any variable to partial out
   ex <- ps$ex
   form <- TRUE
} else {form <- FALSE}
           #  data <- char2numeric(data)   #move to later (01/05/19)
        
    if(is.numeric(y )) y <- colnames(data)[y]
    if(is.numeric(x )) x <- colnames(data)[x]
    if(is.numeric(z )) z <- colnames(data)[z]
  

#check for bad input  
if(any( !(c(y,x,z,ex) %in% colnames(data)) )) {
  cat("\nOops! Variable names are incorrect. Offending items are ", c(y, x, z, ex)[which(!(c(y, x, z, ex) %in% colnames(data)))],"\n")
 stop("I am stopping because the variable names are incorrect.  See above.")}
 
 
 #we might have had formula input thus, make C the covariance matrix
# if(form) C <- cov(data[c(y,x,z,ex)],use="pairwise")
 if(is.null(n.obs)) n.obs <- 0
   numx <- length(x)   #this is the number of predictors (but we should adjust by the number of covariates)   
   numz <- length(z)
   numy <- length(y)
   
   
   #df <- n.obs -1 - numx - length(z) - length(m)    #but this does not take into account the mediating variables
   #note that the x variable includes the intercept and thus uses up one extra df
   df <- n.obs  - numx -numz #We have partialed out z, should we use the df from it?  This is changed 11/26/19 to reduce df for z  
   Cr <- cov2cor(C)  
        	if(!is.null(z)){numz <- length(z)      #partial out the z variables
     	                zm <- C[z,z,drop=FALSE]
     	                za <- C[x,z,drop=FALSE]
     	                zb <- C[y,z,drop=FALSE]
     	                zmi <- solve(zm)
     	                 x.matrix <- C[x,x,drop=FALSE] - za %*% zmi %*% t(za)
     	               if(!part)  y.matrix <- C[y,y,drop=FALSE] - zb %*% zmi %*% t(zb)  #part versus partial (default)
     	                xy.matrix <- C[x,y,drop=FALSE] - za  %*% zmi %*% t(zb)
     	                 C <- cbind(rbind(y.matrix,xy.matrix),rbind(t(xy.matrix),x.matrix))
     	                
     	                 }
   
    if(numx==1) { beta <- solve(C[x,x,drop=FALSE],(C[x,y,drop=FALSE])) 
                 colnames(beta) <- y
        } else {
        beta <- solve(C[x,x],(C[x,y])) }    #this is the same as setCor and is a x * x matrix
        
    if(!is.matrix(beta)) {beta <- matrix(beta,nrow=length(beta))}   #beta is a matrix of beta weights 
    if(is.character(x)) {rownames(beta) <- x}  else {rownames(beta) <- colnames(C)[x]}
    if(is.character(y)) { colnames(beta) <- y} else { colnames(beta) <- colnames(C)[y]}
      
      x.inv <- solve(C[x,x]) #solve x.matrix    #taken from setCor
      yhat <- t(C[x,y,drop=FALSE]) %*% x.inv %*% C[x,y,drop=FALSE]
      resid <- C[y,y]- yhat
      

     if(!std ) {
        df <- n.obs - numx - numz

         Residual.se <- sqrt(diag(resid /df))  #this is the df  n.obs - length(x))
         se <- MSE <- diag(resid )/(df) 
             
               if(length(y) > 1) {SST <- diag(C[y,y] - means[y]^2 * n.obs)} else {SST <- ( C [y,y] - means[y]^2 * n.obs)}
     	       R2 <- (SST - diag(resid)) /SST 
     	       se.beta <- list()
     	        for (i in 1:length(y)) {
     	        se.beta[[i]] <- sqrt(MSE[i] * diag(x.inv))
                 }
                se <- matrix(unlist(se.beta),ncol=numy)
      if(length(y) > 1) {SST <- diag(C [y,y] - means[y]^2 * n.obs)} else {SST <- ( C [y,y] - means[y]^2 * n.obs)}
	       R2 <- (SST - diag(resid)) /SST 
	       
	   
      } else {      R2 <- colSums(beta * C[x,y])/diag(C[y,y,drop=FALSE])   #the standardized case
        uniq <- 1-(1-1/diag(solve(Cr[x,x,drop=FALSE])))  #1- smc
     

        if(n.obs > 2) { # se <- (sqrt((1-R2)/(n.obs-1 - numx-numz)) %*% t(sqrt(1/uniq)))  #these are the standardized se
                        se <- (sqrt((1-R2)/(df)) %*% t(sqrt(1/uniq)))   #setCor uses df = n.obs - numx - 1
                        se <- t( se * sqrt(diag(C[y,y,drop=FALSE])) %*% t(sqrt(1/diag(C[x,x,drop=FALSE]))) )  #But does this work in the general case?

        
                    colnames(se) <- colnames(beta) } else {se <- NA}
                   if(raw) {   #used to compare models  -- we need to adjust this for dfs
                           Residual.se <-   sqrt((1-R2)* df/(df-1)) } else {  #this is a kludge and is necessary to treat the SSR correctly
                       Residual.se <- sqrt((1-R2)/df * (n.obs-1))}
                   }
                   
                
                    if(!any(is.na(se))) { tvalue <- beta/se
                                        # prob <- 2*(1- pt(abs(tvalue),df))
                                          prob <- -2 *  expm1(pt(abs(tvalue),df,log.p=TRUE))
                                         } else {tvalue <- prob <- df <- NA}
  result <- list(beta=beta,se=se, t=tvalue,df=df,prob=prob,R2=R2,SE.resid=Residual.se)
  return(result)       }  		 
#######


