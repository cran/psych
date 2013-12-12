#November 30, 2013
#parts adapted from combn
"splitHalf"<- 
function(r,n.sample=10000,raw=FALSE,covar=FALSE) {
cl <- match.call()
split <- function(o,n) {
A <- B <-  rep(0,n)
A [o] <- B[-o]  <- 1
A[-o] <- B[o] <- 0
AB <- cbind(A,B)
R <- t(AB) %*% r %*% AB
Rab <- R[1,2]/sqrt(R[1,1]*R[2,2])
rab <- 2*Rab/(1+Rab)
result <- list(rab=rab,AB=AB)}

maxrb <- -9999
minrb <- 2
n <- ncol(r)
n2 <- trunc(n/2)
n.obs <- nrow(r)
if(n.obs > n) r <- cov(r,use="pairwise")
if(!covar) r <- cov2cor(r)

e <- 0
m <- n2
h <- m
alpha <- ((sum(r) - tr(r))/sum(r)) * n/(n-1)
sumr <- 0 
x <- seq_len(n)
a <- seq_len(m)

if(is.null(n.sample)  ) {
count <- as.integer(round(choose(n, m)))/2
result <- rep(NA,count)
 #first do the original order
          o <- a
          sp <- split(o,n)
          result[1] <- sp$rab
          maxrb <- sp$rab
         maxAB <- sp$AB
          minrb <- sp$rab
           minAB <- sp$AB
i <- 2 #now, do the rest

 while (i < (count+1)) {
 if (e < n - h) {
                h <- 1L
                e <- a[m]
                j <- 1L
            }
            else {
                e <- a[m - h]
                h <- h + 1L
                j <- 1L:h
            }
             a[m - h + j] <- e + j 
            o <-  x[a]
            sp <- split(o,n)
          result[i] <- sp$rab
          sumr <- sumr+ sp$rab
          if(sp$rab > maxrb) {maxrb  <- sp$rab
                              maxAB <- sp$AB}
        if(sp$rab < minrb) {minrb <- sp$rab
                              minAB <- sp$AB}
i <- i + 1L
}} else {

result <- rep(NA,n.sample)
 for (i in 1:n.sample) {
 o <- sample(n,n2)
 sp <- split(o,n)
result[i] <- sp$rab
          sumr <- sumr+ sp$rab
          if(sp$rab > maxrb) {maxrb <- sp$rab
                              maxAB <- sp$AB}
        if(sp$rab < minrb) { minrb <- sp$rab
                              minAB <- sp$AB}
    }
 }
if(is.null(n.sample)) {meanr <- sumr/count } else {meanr <- sumr/n.sample }
meansp <- 2 * meanr/(1+meanr)
if(raw) {results <- list(maxrb=maxrb,minrb=minrb,maxAB=maxAB,minAB=minAB,meanr=meanr,alpha=alpha,raw = result,Call = cl)
} else {results <- list(maxrb=maxrb,minrb=minrb,maxAB=maxAB,minAB=minAB,meanr=meanr,alpha=alpha,Call=cl)}
class(results) <- c("psych","split")
results
}
