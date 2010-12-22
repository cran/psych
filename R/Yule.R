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
    function(Q,m) {#find the cells for a particular Yule value with fixed marginals
     R1 <- m[1]
     R2 <- m[2]
     C1 <- m[3]
     C2 <- m[4]
     f <- function(x) {(Q- Yule(c(x,R1-x,C1-x,R2-C1+x)))^2}
     xval <- optimize(f,c(0,R1))
     R <- matrix(ncol=2,nrow=2)
     x <- xval$minimum
     R[1,1] <- x
     R[1,2] <- R1-x
     R[2,1] <- C1 -x
     R[2,2] <- R2-C1 + x
     return(round(R))
     }
     
     
     
"Yule2phi" <- 
   function(Q,m) { #find the phi equivalent to a Yule Q with fixed marginals
    t <- Yule.inv(Q,m)
   r.sum <- rowSums(t)
   c.sum <- colSums(t)
   total <- sum(r.sum)
   r.sum <- r.sum/total
   c.sum <- c.sum/total
   v <- prod(r.sum, c.sum)
   phi <- (t[1,1]/total - c.sum[1]*r.sum[1]) /sqrt(v)
 return(phi)}
 
 
"Yule2poly" <-
    function(Q,m) { #find the phi equivalent to a Yule Q with fixed marginals
    
    t <- Yule.inv(Q,m)
    r <- tetrachoric(t)
    return(r)
    }
    
    
    
    

