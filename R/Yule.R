"Yule" <- 
 function(x,Y=FALSE) {
  # if(!is.matrix(x)) {stop("x must be a matrix")}
   stopifnot(prod(dim(x)) == 4 || length(x) == 4)
    if (is.vector(x)) {
     x <- matrix(x, 2)}
   a  <- x[1,1]
   b <- x[1,2]
   c <- x[2,1]
   d <- x[2,2]
  if (Y) {Yule <- (sqrt(a*d) - sqrt(b*c))/(sqrt(a*d)+sqrt(b*c))} else {Yule <- (a*d- b*c)/(a*d + b*c)}
   return(Yule) #Yule Q
   }
   
 
"Yule.inv" <- 
    function(Q,m,n=NULL) {#find the cells for a particular Yule value with fixed marginals
     if (length(m) > 2) { #this old way is for one correlation at a time with all 4 marginals - 
     R1 <- m[1]
     R2 <- m[2]
     C1 <- m[3]
     C2 <- m[4] } else {   #the better way is to specify just two marginals- allows better interface with matrices
       if(is.null(n)) {n<- 1}  #do it as percentages
        
        R1 <- m[1] 
        R2 <- n -R1
        C1 <- m[2]
        C2 <- n - C1}
     
     f <- function(x) {(Q- Yule(c(x,R1-x,C1-x,R2-C1+x)))^2}
     xval <- optimize(f,c(min(R1-min(C2,R1),C1-min(R2,C1)),min(R1,C1)))
     R <- matrix(ncol=2,nrow=2)
     x <- xval$minimum
     R[1,1] <- x
     R[1,2] <- R1-x
     R[2,1] <- C1 -x
     R[2,2] <- R2-C1 + x
     return(R)
     }
     
"Yule2poly" <-
    function(Q,m,n=NULL,correct=TRUE) { #find the phi equivalent to a Yule Q with fixed marginals
    .Deprecated("Yule2tet",msg="Yule2poly has been replaced by Yule2tet, please try again") 
    t <- Yule.inv(Q,m,n=n)
    r <- tetrachoric(t,correct=correct)$rho
    return(r)
    }

"Yule2tet" <-
    function(Q,m,n=NULL,correct=TRUE) { #find the phi equivalent to a Yule Q with fixed marginals
    t <- Yule.inv(Q,m,n=n)
    r <- tetrachoric(t,correct=correct)$rho
    return(r)
    }
    
"Yule2phi" <- 
   function(Q,m,n=NULL) { #find the phi equivalent to a Yule Q with fixed marginals
    t <- Yule.inv(Q,m,n=n)
  return(phi(t,digits=15))}
  
  
  
"Yule2tetra" <-
function(Q,m,n=NULL,correct=TRUE) {
if(!is.matrix(Q) && !is.data.frame(Q)) {result <- Yule2tet(Q,c(m[1],m[2]),n=n,correct=correct) } else {
nvar <- nrow(Q)
if(nvar !=ncol(Q)) {stop('Matrix must be square')}
if (length(m) !=nvar) {stop("length of m must match the number of variables")}
result <- Q
for(i in 2:nvar) {
  for (j in 1:(i-1)) { 
  result[i,j] <- result[j,i] <- Yule2tet(Q[i,j],c(m[i],m[j]),n,correct=correct)
  }
}}
return(result) }

"YuleBonett" <- 
 function(x,c=1,bonett=FALSE,alpha=.05) {
  # if(!is.matrix(x)) {stop("x must be a matrix")}
   stopifnot(prod(dim(x)) == 4 || length(x) == 4)
    if (is.vector(x)) {
     x <- matrix(x, 2)}
     p <- x/sum(x)
    C <- c
   a  <- p[1,1]
   b <- p[1,2]
   c <- p[2,1]
   d <- p[2,2]
   Rs <- rowSums(p)
   Cs <- colSums(p)
   if(bonett) {C <- .5 - (.5 - min(Rs,Cs))^2 }
   ad <- (a*d)^C
   bc <- (b*c)^C
   Yule <- (ad-bc)/(ad+bc) 
   Ystar <- 
   #See Bonett 2007 p 434
   w <- (x[1,1] + .1)* (x[2,2] + .1)/((x[2,1] + .1)* (x[1,2] + .1)) #OR estimate
    Ystar <- (w^C -1)/(w^C+1)   #this is the small cell size adjusted value
   	vlQ <- (1/(x[1,1] + .1) + 1/ (x[2,2] + .1) + 1/(x[2,1] + .1)+ 1/ (x[1,2] + .1))
	vQ <-  (C^2/4)*(1-Ystar^2)^2 *vlQ   #equation 9  
  tanhinv <- atanh(Ystar)
  upper <- tanh(atanh(Ystar) + qnorm(1-alpha/2) *sqrt( (C^2/4) *vlQ))
   lower <- tanh(atanh(Ystar) - qnorm(1-alpha/2) *sqrt( (C^2/4) *vlQ))
   result <- list(rho=Ystar,se=sqrt(vQ), upper=upper, lower=lower)
   class(result) <- c("psych","Yule")
   return(result) #Yule Q, Y, or generalized Bonett
   }
   
   
   
#Find a correlation matrix of Yule correlations
"YuleCor" <- function(x,c=1,bonett=FALSE,alpha=.05) {
cl <- match.call()
nvar <- ncol(x)
rho <- matrix(NA,nvar,nvar)
ci <-  matrix(NA,nvar,nvar)
zp <- matrix(NA,nvar,nvar)
for(i in 1:nvar) {
  for(j in 1:i) {
    YB <- YuleBonett(table(x[,i],x[,j]),c=c,bonett=bonett,alpha=alpha)
    rho[i,j] <- rho[j,i] <- YB$rho
    ci[i,j] <- YB$lower
    ci[j,i] <- YB$upper
    zp[i,j] <- YB$rho/YB$se
    zp[j,i] <- 1-pnorm(zp[i,j])
    
    }}
colnames(rho) <- rownames(rho) <- colnames(ci) <- rownames(ci) <- colnames(x)
result <- list(rho=rho,ci=ci,zp=zp,Call=cl)
class(result) <- c("psych","yule")
return(result)
}

  
    
    

