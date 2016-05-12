"correct.cor" <-
function(x,y) { n=dim(x)[1]   
        { diag(x) <- y 
        if (n> 1)  {
        for (i in 2:n) {
           k=i-1
           for (j in 1:k) {
              x[j,i] <- x[j,i]/sqrt(y[i]*y[j])  }   #fix the upper triangular part of the matrix
             }}
           return(x)  }}
           
"rangeCorrection" <- function(r,sdu,sdr,sdxu=NULL,sdxr=NULL,case=2) {
if (!is.null(sdxu)) case <- 4  #
switch(case,
{ result <-  sqrt(1-(sdr^2/sdu^2) *(1-r^2))},
{ result <- (  r * sdu/(sdr* sqrt(1-r^2 + r^2*(sdu^2/sdr^2))))},
{result <- NULL},
{result <- r * (sdr/sdu)*(sdxr/sdxu) + sqrt((1-(sdr/sdu)^2) * (1- (sdxr/sdxu)^2 ))  }
)
return(result)
}



#Find the Kaiser - Meyer -Olkin criterion
#note that the correct formula is in Kaiser 1974, not 1970
"KMO" <-
function(r) {
cl <- match.call()
if(nrow(r) > ncol(r)) r <- cor(r,use="pairwise")
  Q <- try(solve(r))
  if(class(Q) == as.character("try-error")) {message("matrix is not invertible, image not found")
        Q <- r}
 S2  <- diag(1/diag(Q))
 IC <- S2 %*% Q %*% S2
 Q <- Image <-  cov2cor(Q) 
 diag(Q) <- 0
 diag(r) <- 0
 sumQ2 <- sum(Q^2)
 sumr2 <- sum(r^2)
 MSA <- sumr2/(sumr2 + sumQ2)
 MSAi <- colSums(r^2)/(colSums(r^2) + colSums(Q^2))
 results <- list(MSA =MSA,MSAi = MSAi, Image=Image,ImCov = IC,Call=cl)
 class(results) <- c("psych","KMO")
 return(results)
 }