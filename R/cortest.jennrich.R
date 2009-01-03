"cortest.jennrich" <- 
 function(R1,R2,n1=NULL, n2=NULL) {
 p <- dim(R1)[2]
if(dim(R1)[1] != p) { n1 <- dim(R1)[1]
                     R1 <- cor(R1,use="pairwise")
                    warning ("R1 matrix was not square, correlations found") 
                    }
 if(dim(R2)[1] != dim(R2)[2] ) {n2 <- dim(R2)[1]
                     R2 <- cor(R2,use="pairwise") 
                    warning ("R2 matrix was not square, correlations found") }
if(!is.matrix(R1) ) R1 <- as.matrix(R1)  #converts data.frames to matrices if needed
if(!is.matrix(R2) ) R2 <- as.matrix(R2)

if (dim(R1)[2] != dim(R2)[2]) stop("correlation matrices M and S must be of the same size!")
 if(is.null(n2)) n2 <- n1
 if(!is.null(n1) & !is.null(n2)) c <- n1*n2/(n1+n2) else c <- 1
 
 R <- (n1*R1+n2*R2)/(n1+n2)   #matrix of average values

 S <- R * R    #squared values of averaged correlations
 S.inv <- solve(S)
 R.inv <- solve(R)
 R.diff <- R1 - R2
 Z <- sqrt(c) * R.inv %*% R.diff
 chi2 <- tr(Z%*%t(Z))/2 - t(diag(Z)) %*% S.inv %*% diag(Z)
 chi2 <- chi2[1,1]
 p <- dim(R1)[1]
 df <- p*(p-1)/2 
 results <- list(chi2 =chi2,prob =pchisq(chi2,df,lower.tail=FALSE))
 return(results)
 }
 #Jennrich (1970) p 908
 