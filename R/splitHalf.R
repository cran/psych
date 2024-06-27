#November 30, 2013
#parts adapted from combn
#modified 6/20/21 to just do the unique splits if < 15,000 *i.e. a 16 variable problem
#modified 9/7/21 to handle keys as lists of variable names with negative signs
"splitHalf"<- 
function(r,raw=FALSE,brute=FALSE,n.sample=15000,covar=FALSE,check.keys=TRUE,key=NULL,ci=.05,use="pairwise") {
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

#if(!is.null(key)) {r <- r[,selectFromKeys(key)]}  #made this again 12/23/23
if(!is.null(colnames(r))) {v.names <- colnames(r)} else {v.names <- colnames(r)<- paste0("V",1:NCOL(r))}

keys <- key   
maxrb <- -9999
minrb <- 2
n <- ncol(r)
n2 <- trunc(n/2)
n.obs <- nrow(r)

#if(n.obs > n) { r <- cov(r,use=use)}    # fails if n.obs < n
if(!isCovariance(r))  { r <- cov(r,use=use)}   
 if(!covar) r <- cov2cor(r) 
 
if(check.keys && is.null(keys)) {
            p1 <- principal(r,covar=covar)
            if(any(p1$loadings < 0)) warning("Some items were negatively correlated with total scale and were automatically reversed.")
            keys <- 1- 2* (p1$loadings < 0 ) 
            }  #keys is now a vector of 1s and -1s
            
         if(is.null(keys)) {keys <- rep(1,n)} else {     #in case keys is a list, we need to tranform to a keys vector 
            if(is.character(keys))  {temp <- rep(1,n)
                  temp[which(strtrim(keys,1)=="-")] <- -1
                  keys <- temp } 
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
                     colnames(r) <- paste(signkey,colnames(r),sep="")

e <- 0
m <- n2
h <- m
sumr <- sum(r)
 off <- r
diag(off) <- 0
sum.off <- sum(off)
sumsq.off <- sum(off^2)
sum.off <- sumr - tr(r)
alpha <- (sum.off/sumr) * n/(n-1)
tsmc <- sum(smc(r,covar=covar))
lambda6 <- (sum.off + tsmc)/sumr
lambda2 <- (sum.off+ sqrt(sumsq.off*n/(n-1)))/sumr
result <- NULL
med.r <- median(r[lower.tri(r)],na.rm=TRUE)  #find the median correlation
av.r <- mean(r[lower.tri(r)],na.rm=TRUE) 
sumr <- 0 
x <- seq_len(n)
a <- seq_len(m)
#count <- as.integer(round(choose(n, m)))/2 #doesn't work for very large values of n.
#count <- round(choose(n, m))/2
count <- round(choose(n, m))  #do all splits because there is a problem with odd values    10/21/23
if(brute || ((count <= n.sample))) {    #brute force -- try all combinations
brute <- TRUE 
 if(raw) result <- rep(NA,count)  #keep the results if raw
 #first do the original order
        o <- a
        sp <- split(o,n)
        if(raw) result[1] <- abs(sp$rab)
        sp$rab <- abs(sp$rab)
        maxrb <- sp$rab
        maxAB <- sp$AB
        minrb <- sp$rab
        minAB <- sp$AB
        sumr <- sp$rab   #added to get the first one 10/6/18 in response to bug report by Wes Bonifay
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
         if(raw) result[i] <- abs(sp$rab )
           sp$rab <- abs(sp$rab)  #conrolling for a rare case
            sumr <- sumr+ sp$rab  # summing all of the split halfs to eventually report mean splithalf
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
  #  )    #if using mclapply
 }
 #now 
if(brute) {meanr <- sumr/count } else {meanr <- sumr/n.sample }


kd <- diag(key.d)    #reverse them so we can use the keys
maxAB =maxAB * kd
minAB = minAB * kd
rownames(maxAB) <- rownames(minAB) <- v.names
maxAB <- keys2list(maxAB)
minAB <- keys2list(minAB)
if(!anyNA(result))  {
ci <- quantile(result,c(ci/2,.5, 1 - ci/2))} else {ci <- rep(NA,3) }
if(raw) {
results <- list(maxrb=maxrb,minrb=minrb,maxAB=maxAB,minAB=minAB,meanr=meanr,av.r=av.r,med.r=med.r, alpha=alpha,lambda2 = lambda2, lambda6=lambda6,raw = result,ci=ci,covar=covar,Call = cl)
} else {results <- list(maxrb=maxrb,minrb=minrb,maxAB=maxAB,minAB=minAB,meanr=meanr,av.r=av.r,med.r=med.r,alpha=alpha,lambda2 = lambda2,lambda6=lambda6,ci=ci,covar=covar, Call=cl)}
class(results) <- c("psych","split")
return(results)
}
