matrix.na <- function (x, y) {
nvar <- ncol(x)
n.obs = nrow(x)
ny <- ncol(y)
if(nvar!= NROW(y)) stop
result <- matrix(NA,ncol=ny,nrow=n.obs)
for (k in 1:n.obs) {
  for (i in 1:ny ) {
   result[k,i] <- mean(x[k,] * y[,i],na.rm=TRUE)}
   }
result
}

#an attempt to do fast matrixlike operation with missing data\

matrix.na.z		 <- function(x,y,scale=TRUE) {
#first get tranpose x to make multiplication work

nvar <- ncol(x)
if(nvar != nrow(y) ) stop("matrices are not compatible")#matrices are not compatible

if(scale) x <- scale(x) #zero center and standaridize
tx <- t(x) #we want to do  it on the transposed matrix
ny <- ncol(y)
result <- matrix(NA,nrow = nrow(x),ncol= ncol(y))
result <- apply(y,2,function(x ) colSums(x * tx,na.rm=TRUE))
return((result))

}

"impute.na" <- function(x,impute="mean") {
  miss <- which(is.na(x),arr.ind=TRUE)
   if(impute=="mean") {
       		item.means <- colMeans(x,na.rm=TRUE)   #replace missing values with means
       		x[miss]<- item.means[miss[,2]]} else { 
       		item.med   <- apply(items,2,median,na.rm=TRUE) #replace missing with medians
        	x[miss]<- item.med[miss[,2]]} 
    return(x)
 }
 
 
  
score.na <- function ( keys,r,cor=TRUE,smooth=FALSE) {#score a matrix with missing correlations
covar <- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))*x,na.rm=TRUE))
count <- t(keys) %*% keys  #counts the number of items/key
bad <- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(is.na(r)*abs(x),na.rm=TRUE))*abs(x),na.rm=TRUE)) #the number missin
good <- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(!is.na(r)*abs(x),na.rm=TRUE))*abs(x),na.rm=TRUE)) #the number not missing

result <- (covar - count)/(good-count) *(good-count + bad)+ count
if(cor) result <- cov2cor(result)
if(smooth) result<- cor.smooth(result)
return(result) #the covariance matrix  use cov2cor
} 