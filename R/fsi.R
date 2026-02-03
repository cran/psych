 #factor score indeterminancy  following Nicewander  2022
 #and Grice 2001
fsi <- function(f,phi=NULL,r=NULL,Grice=FALSE,short=TRUE) {
if(Grice){R2 <- gricef(f=f,phi=phi,r=r)
     Hn <- NULL} else {
if(!is.matrix(f)) f <- as.matrix(f)
H2 <- rowSums(f^2)
H2[H2>1] <- 1
fe <- diag(sqrt(1-H2))
L <- cbind(f,fe)
if(is.null(phi)) {H <- Pinv(L) %*% L } else {
 fr <-  f %*% phi
 H2 <- rowSums(f^2)
 H2[H2>1] <- 1
fe <- diag(sqrt(1-H2))
L <- cbind(fr,fe)
 H <- Pinv(L) %*% L}
 R2 <- diag(H) 
 Hn <- H[1:ncol(f),1:ncol(f)]}
 names(R2) <- colnames(f) 
if(short) {result <- R2[1:ncol(f)]} else {result <- list(R2=R2, H = Hn)}
return(result)   #these are the R2 of the  factors with the scores
#note that fungible:fsIndeterminacy  returns r
}  

gricef <- function(f,phi=NULL,r=NULL) {
    if(is.null(phi) ){phi <- diag(ncol(f))}
      S <- f %*% phi 
      w  <- Pinv(r) %*% S 
      R2 <- diag(t(w) %*% S)
      R2[R2 >1] <- 1 
      return(R2)
      }