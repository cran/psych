#November 30, 2013
#parts adapted from combn
"splitHalf"<- 
function(r,raw=FALSE,brute=FALSE,n.sample=10000,covar=FALSE,check.keys=TRUE,key=NULL) {
cl <- match.call()
split <- function(o,n) {
A <- B <-  rep(0,n)
A [o] <- B[-o]  <- 1
A[-o] <- B[o] <- 0
AB <- cbind(A,B)
R <- t(AB) %*% r %*% AB
Rab <- R[1,2]/sqrt(R[1,1]*R[2,2])
#rab <- 2*Rab/(1+Rab)
rab <- 4*R[1,2]/sum(R)
result <- list(rab=rab,AB=AB)}

keys <- key   
maxrb <- -9999
minrb <- 2
n <- ncol(r)
n2 <- trunc(n/2)
n.obs <- nrow(r)
if(n.obs > n) { r <- cov(r,use="pairwise")}
 if(!covar) r <- cov2cor(r) 
 
if(check.keys && is.null(keys)) {
            p1 <- principal(r)
            if(any(p1$loadings < 0)) warning("Some items were negatively correlated with total scale and were automatically reversed.")
            keys <- 1- 2* (p1$loadings < 0 ) 
            }  #keys is now a vector of 1s and -1s
            
         if(is.null(keys)) {keys <- rep(1,n)} else {  
         			keys<- as.vector(keys) 
         			if(length(keys) < n) {temp <- keys  #this is the option of keying just the reversals
         			                         keys <- rep(1,n)
         			                         names(keys) <- colnames(r)
         			                         keys[temp] <- -1}
         			                   } 
        			 key.d <- diag(keys)
                     r <- key.d %*% r %*% key.d
                     signkey <- strtrim(keys,1)
            		 signkey[signkey=="1"] <- ""
                     colnames(r) <- paste(colnames(r),signkey,sep="")

e <- 0
m <- n2
h <- m
sumr <- sum(r)
alpha <- ((sumr - tr(r))/sumr) * n/(n-1)
tsmc <- sum(smc(r))
lambda6 <- (sumr - tr(r) + tsmc)/sumr


sumr <- 0 
x <- seq_len(n)
a <- seq_len(m)
#count <- as.integer(round(choose(n, m)))/2 #doesn't work for very large values of n.
count <- round(choose(n, m))/2

if(brute || ((count <= n.sample) && !raw)) {    #brute force -- try all combinations
brute <- TRUE 
 if(raw) result <- rep(NA,count)  #keep the results if raw
 #first do the original order
        o <- a
        sp <- split(o,n)
        if(raw) result[1] <- sp$rab
        maxrb <- sp$rab
        maxAB <- sp$AB
        minrb <- sp$rab
        minAB <- sp$AB
i <- 2 #now, do the rest
 while (i < (count+1)) { #adapted (taken) from combn
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
         if(raw) result[i] <- sp$rab
          sumr <- sumr+ sp$rab
          if(sp$rab > maxrb) {maxrb  <- sp$rab
                              maxAB <- sp$AB}
        if(sp$rab < minrb) {minrb <- sp$rab
                              minAB <- sp$AB}
i <- i + 1L
}} else {   #sample the alternatives

result <- rep(NA,n.sample)
sumr <- 0
for (i in 1:n.sample) {
#result <- mclapply(1:n.sample,{   #use mclapply to allow for parallelism
 	o <- sample(n,n2)
 	sp <- split(o,n)
	if(raw) result[i] <- sp$rab
	#if(raw) result <- sp$rab   #if mclapply
    sumr <- sumr+ sp$rab
    if(sp$rab > maxrb) {maxrb <- sp$rab
                              maxAB <- sp$AB}
    if(sp$rab < minrb) { minrb <- sp$rab
                              minAB <- sp$AB}
                        
    }
    #)    #if using mclapply
 }
 #now 
if(brute) {meanr <- sumr/count } else {meanr <- sumr/n.sample }
meansp <- 2 * meanr/(1+meanr)

kd <- diag(key.d)    #reverse them so we can use the keys
maxAB = maxAB * kd
minAB = minAB * kd

if(raw) {results <- list(maxrb=maxrb,minrb=minrb,maxAB=maxAB,minAB=minAB,meanr=meanr,alpha=alpha,lambda6=lambda6,raw = result,Call = cl)
} else {results <- list(maxrb=maxrb,minrb=minrb,maxAB=maxAB,minAB=minAB,meanr=meanr,alpha=alpha,lambda6=lambda6,Call=cl)}
class(results) <- c("psych","split")
return(results)
}
