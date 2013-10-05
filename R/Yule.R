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

    
    
    

