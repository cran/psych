"cortest.mat" <- 
function(R1,R2=NULL,n1=NULL,n2 = NULL) {
cl <- match.call()
p <- dim(R1)[2]
if(dim(R1)[1] != p) { n1 <- dim(R1)[1]
                    R1 <- cor(R1,use="pairwise")
                    warning ("R1 matrix was not square, correlations found") 
                    }
 if(!is.matrix(R1)) R1 <- as.matrix(R1)  #in case R1 is a data.frame
 
 if (is.null(n1)) {warning("n1 not specified, set to 100")
            n1 <-100}
 if (is.null(R2)) {message("Bartlett's test of is R = I") 
                  detR1 <- det(R1)
                chi2  <- -log(detR1) *(n1 -1 - (2*p + 5)/6)
                  df <- p * (p-1)/2
                  pval <- pchisq(chi2,df,lower.tail=FALSE)
                  n.obs <- n1
                } else {
                       if(dim(R2)[1] != dim(R2)[2] ) {n2 <- dim(R2)[1]
                         R2 <- cor(R2,use="pairwise") 
                    warning ("R2 matrix was not square, correlations found") }
                    
if (p != dim(R2)[2]) stop("correlation matrices R1 and R2 must be of the same size!")
 if(!is.matrix(R2)) R2 <- as.matrix(R2)  #in case R1 is a data.frame
R1.inv <- solve(R1)    #inverse of R
R2.inv <- solve(R2)    #inverse of R

R.inv.2 <- R1.inv %*% R2 #inverse of R1 times R2
R.inv.1 <- R2.inv %*% R1 #inverse of R2 times R1

E1 <- .5*(sum((diag(R.inv.2))) -log(det(R.inv.2)) - p ) #likelihood

E2 <- .5*(sum((diag(R.inv.1))) -log(det(R.inv.1)) - p ) #likelihood

df1 <- p * (p-1)/2
df <- 2*df1
 
if (is.null(n2)) {n2 <- n1}
n.obs <- min(n1,n2)
chi21 <- E1 * (n1-1-(2*p - 5)/6)
chi22 <- E2 * (n2-1-(2*p - 5)/6)
chi2 <- chi21 + chi22} 
 results <- list(chi2 =chi2,prob =pchisq(chi2,df,lower.tail=FALSE), df= df,n.obs=n.obs,Call =cl)
 class(results) <- c("psych","cortest")
 return(results)
}
