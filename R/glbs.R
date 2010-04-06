"glb.fa" <- 
function(r,key= NULL){

 cl <- match.call()   #for eventual fancy printing with the call listed
 nvar <- dim(r)[2]    #find a correlation matrix if using a data matrix
if(dim(r)[1] != dim(r)[2]) {r <- cor(r,use="pairwise")}  else {
if(!is.matrix(r)) r <- as.matrix(r)
r <- cov2cor(r)}  #make sure it is a correlation matrix not a covariance or data matrix
      if(is.null(colnames(r))) {  rownames(r) <- colnames(r) <- paste("V",1:nvar,sep="") }
      if (!is.null(key)) { key <- as.vector(key)
                          r <- diag(key) %*% r %*% diag(key)
                       
                           flip <- FALSE   #we do this if we specify the key
                           } else {key <- rep(1,nvar) }
   nv <- dim(r)[1]         #how many variables 
   f1 <- fa(r) #factor it  #first find the eigen values of the factor model
   nf <-length(which(f1$values > 0))  #how many are real 
   df <- nv * (nv-1)/2  - nf*nv + nf*(nf-1)/2  #check for degrees of freedom   
   if (df <0 ) nf <- nf-1
   fn <- fa(r,nf,rotate="none")
   rr <- r
   diag(rr) <- fn$communality
   glb <- sum(rr)/sum(r)
   return(list(glb=glb,communality = fn$communality,numf = nf,Call=cl))
   }
    