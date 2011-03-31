
"sdt" <- 
 function(x) {
  # if(!is.matrix(x)) {stop("x must be a matrix")}
   stopifnot(prod(dim(x)) == 4 || length(x) == 4)
    if (is.vector(x)) {
     x <- matrix(x, 2)}
   a  <- x[1,1]
   b <- x[1,2]
   c <- x[2,1]
   d <- x[2,2]
   colnames(x) <- c("signal","noise")
   rownames(x) <- c("Signal","Noise")
   
 T <- sum(x)
 signals <- a+b
 H <- a/signals
 noise <- c+d
 F <- c/noise
 dprime <- qnorm(H) - qnorm(F)
 beta <- exp(((qnorm(F)^2 - qnorm(H))^2/2)) 
 Beta <- dnorm(qnorm(H))/dnorm(qnorm(F))
 Aprime <- .5 + sign((H-F) *((H-F)^2 + abs(H-F))/(4*max(H,F)-4 * H*F))
 aprime <- .5 +(H-F)*(1-H-F)/(4*H*(1-F))
  C <- -(qnorm(H)+qnorm(F))/2 
  Bprime <- sign(H-F) * ((H*(1-H) - F*(1-F))/((H*(1-H) + F*(1-F)) ))
 sdt <- list(x=x,dprime=dprime,Beta=Beta,lnBeta=log(Beta),beta=beta,lnbeta =log(beta),Aprime=Aprime,C=C,Bprime=Bprime,aprime=aprime)
return(sdt)
 }
 
 
 beta.w <- function(hit,fa) {
 zhr <- qnorm(hit)
 zfar <- qnorm(fa) 
 exp(-zhr*zhr/2+zfar*zfar/2)
}

aprime <-function(hit,fa) { a<-1/2+((hit-fa)*(1+hit-fa) / (4*hit*(1-fa))) 
b<-1/2-((fa-hit)*(1+fa-hit) / (4*fa*(1-hit)))
a[fa>hit] <-b[fa>hit]
a[fa==hit]<-.5 
a
}
 
 test.sdt <- function(n) {
 hits <- 1:n
 misses <- (n-1):0
 fa <- sample(n,n,replace=TRUE)
# fa <- n/4
 cr <- n-fa
 test.sdt <- data.frame(hits,misses,fa,cr)
 }
 
 
 "sdt" <- 
 function(x) {
 if(is.null(dim(x))) {# the case of a single problem
 if (length(x) == 2) {x <- matrix(x,nrow=1)
     H <- x[1]
     F <- x[2]
    colnames(x) <- c("Hit Rate","False Alarmes")} else {
     x <- matrix(x,nrow=1)
    # H <- x[,1]/(x[,1] + x[,2])
     F <- x[,3]/(x[,3] + x[,4])
    x[(x[,2] < 1)] <- .5
    H <- x[,1]/(x[,1] + x[,2])
    if (x[,3] < 1) {F <- x[,3]/(x[,3] + x[,4] + .5)}
      colnames(x) <- c("signal","noise","Signal","Noise")}} else {
   if(dim(x)[2] ==2) { H <- x[,1]
     F <- x[,2] 
      colnames(x) <- c("Hit Rate","False Alarmes")} else { #the vector case
     x[(x[,2] < 1),2] <- .5
     x[(x[,3] < 1),3] <- .5
     
      H <- x[,1]/(x[,1] + x[,2])
      F <- x[,3]/(x[,3] + x[,4])
   
   colnames(x) <- c("signal","noise","Signal","Noise")}
}
 
 dprime <- qnorm(H) - qnorm(F)
 beta <- dnorm(qnorm(H))/dnorm(qnorm(F)) #definitional
 cprime <-  -.5*(qnorm(H) + qnorm(F))/dprime   #macmillan
 beta.w <- exp(-qnorm(H)^2/2+qnorm(F)/2 ) #pallier  = beta    
 Aprime <- .5 +(H-F)*(1+H-F)/(4*H*(1-F)) #Grier
 Bprime <- sign(H-F) * ((H*(1-H) - F*(1-F))/((H*(1-H) + F*(1-F)) ))
 bprime <- ((1-H)*(1-F)-H*F)/((1-H)*(1-F)+H*F)
 sdt <- data.frame(x=x,dprime=dprime,beta=beta,cprime=cprime,lnBeta=log(beta),Aprime=Aprime,Bprime=Bprime,brime=bprime,beta.w = beta.w)
return(sdt)
 }