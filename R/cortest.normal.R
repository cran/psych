"cortest.normal" <- 
function(R1,R2=NULL, n1=NULL,n2=NULL,fisher=TRUE) {
cl <- match.call()
if (dim(R1)[1] != dim(R1)[2]) {n1 <- dim(R1)[1] 
                             message("R1 was not square, finding R from data")
                             R1 <- cor(R1,use="pairwise")}
 
if(!is.matrix(R1) ) R1 <- as.matrix(R1)  #converts data.frames to matrices if needed

p <- dim(R1)[2]
if(is.null(n1)) {n1 <- 100 
                warning("n not specified, 100 used") }
if(is.null(R2)) { if(fisher) {R <- 0.5*log((1+R1)/(1-R1))
                              R <- R*R} else {R <- R1*R1}
                 diag(R) <- 0
                 E <- (sum(R*lower.tri(R)))
                 chisq <- E *(n1-3)
                 df <- p*(p-1)/2
                 p.val <- pchisq(chisq,df,lower.tail=FALSE)
    } else {         #end of 1 matrix test
    if (dim(R2)[1] != dim(R2)[2]) {n2 <- dim(R2)[1] 
                             message("R2 was not square, finding R from data")
                             R2 <- cor(R2,use="pairwise")}
      if(!is.matrix(R2) ) R2 <- as.matrix(R2)

                             
      if(fisher) { 
                  R1 <- 0.5*log((1+R1)/(1-R1)) 
                  R2 <-  0.5*log((1+R2)/(1-R2))
                  diag(R1) <- 0
                  diag(R2) <- 0 }
        R <-  R1 -R2   #direct difference 
        R <- R*R
        if(is.null(n2)) n2 <- n1
        n <- (n1*n2)/(n1+n2)
         E <- (sum(R*lower.tri(R)))
                 chisq <- E *(n-3)
                 df <- p*(p-1)/2
                 p.val <- pchisq(chisq,df,lower.tail=FALSE)
      }
   result <- list(chi2=chisq,prob=p.val,df=df,Call=cl)
   class(result) <- c("psych", "cortest")
   return(result)
    }
#version of 2008


#the following is done for non symmetric matrices with the same logic
#version of August 28,2011
#not yet ready for prime time
"cortest.normal1" <- 
function(R1,R2=NULL, n1=NULL,n2=NULL,fisher=TRUE) {
cl <- match.call()

 
if(!is.matrix(R1) ) R1 <- as.matrix(R1)  #converts data.frames to matrices if needed
if(!is.matrix(R2) ) R2 <- as.matrix(R2)  #converts data.frames to matrices if needed

r <- dim(R1)[1]
c <- dim(R1)[2]
   
                  R1 <- 0.5*log((1+R1)/(1-R1)) 
                  R2 <-  0.5*log((1+R2)/(1-R2))
                 
        R <-  R1 -R2   #direct difference 
        R <- R*R
        if(is.null(n2)) n2 <- n1
        n <- (n1*n2)/(n1+n2)
         E <- sum(R)
                 chisq <- E *(n-3)
                 df <- r*c
                 p.val <- pchisq(chisq,df,lower.tail=FALSE)

   result <- list(chi2=chisq,prob=p.val,df=df,Call=cl)
   class(result) <- c("psych", "cortest")
   return(result)
    }




