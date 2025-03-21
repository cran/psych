#added SRMR and RMSEA  1/11/25
"cortest" <- 
function(R1,R2=NULL, n1=NULL,n2=NULL,fisher=TRUE,cor=TRUE, method = "pearson",use="pairwise") {
cl <- match.call()

if ((NROW(R1)!= NCOL(R1))  & cor) {n1 <- NROW(R1)
                            # message("R1 was not square, finding R from data")
                             R1 <- cor(R1,use=use,method=method)}
 
if(!is.matrix(R1) ) R1 <- as.matrix(R1)  #converts data.frames to matrices if needed

p <- NCOL(R1)
if(is.null(n1)) {n1 <- 100 
                warning("n not specified, 100 used") }
if(is.null(R2)) { if(fisher) {R <- 0.5*log((1+R1)/(1-R1))
                              R2 <- R*R} else {R2 <- R1*R1}
                 if(cor) {diag(R2)  <- 0
                 E <- (sum(R2*lower.tri(R2)))
                 z <- sum(R2*lower.tri(R2))  
                 df <- p*(p-1)/2
                 } else {
                     E <- sum(R2)
                     z <- sum(R1^2) 
                     df <- ncol(R1) * nrow(R1)}
                  chisq <- E *(n1-3)
                  n <- n1
                 z <- z/df
                 
                 p.val <- pchisq(chisq,df,lower.tail=FALSE)
    } else {         #end of 1 matrix test
    if ((dim(R2)[1] != dim(R2)[2]) & cor)  {n2 <- dim(R2)[1] 
                             message("R2 was not square, finding R from data")
                             R2 <- cor(R2,use=use, method=method)}
      if(!is.matrix(R2) ) R2 <- as.matrix(R2)

                             
      if(fisher) { 
                  R1 <- 0.5*log((1+R1)/(1-R1)) 
                  R2 <-  0.5*log((1+R2)/(1-R2)) 
                  if(cor) {diag(R1) <- 0
                  diag(R2) <- 0} }
        R <-  R1 - R2   #direct difference 
        R2 <- R*R
        if(is.null(n2)) n2 <- n1
        n <- 2*(n1*n2)/(n1+n2)  # harmonic sample size 
        if(cor) { E <- (sum(R2*lower.tri(R2)))  #just count the lower diagonal elements
                 chisq <- E *(n-3)
                 df <- p*(p-1)/2
                 z <- sum(R2*lower.tri(R2)) / df} else {E <- sum(R2)
                   chisq <- E * (n-3)
                   df <- ncol(R2) * nrow(R2)
                   z <- sum(R2) / df
                   }
                 p.val <- pchisq(chisq,df,lower.tail=FALSE)
      }
   # if (is.null(n2) ) z <- NULL
   
    SRMR <- sqrt(z)
    RMSEA <- sqrt(max(chisq/(df* n) - 1/(n-1), 0))        #this is x2/(df*N ) -  1/(N-1)   #fixed 4/5/19
   result <- list(chi2=chisq,prob=p.val,df=df,z,RMSEA=RMSEA,SRMR=SRMR,z =z,Call=cl)
   class(result) <- c("psych","cortest")
   return(result)
    }


#SMR amd RMSEA added January 10, 2024
#version of June 25, 2008
#revised October 12, 2011 to allow non-square matrices


test.cortest <- function(R=NULL,n.var=10,n1=100,n2=1000,n.iter=1000) {
if(is.null(R)) R <- diag(1,n.var)
summary <- list()
for(i in 1:n.iter) {
x <- sim.correlation(R,n1)
if(n2 >3 ) {
y <- sim.correlation(R,n2)
summary[[i]] <- cortest(x,y,n1=n1,n2=n2)$prob
} else {summary[[i]] <- cortest(x,n1=n1)$prob }
}
result <- unlist(summary)
return(result)
}

