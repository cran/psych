#corrected May 7, 2007 
"skew" <- 
function (x, na.rm = TRUE) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
        sdx <- sd(x,na.rm=na.rm)
        mx <- mean(x)
        skewer <- sum((x - mx)^3)/(length(x) * sd(x)^3) 
        } else {
    
    skewer <- rep(NA,dim(x)[2])
    if (is.matrix(x)) {mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
    sdx <- sd(x,na.rm=na.rm)
    for (i in 1:dim(x)[2]) {
    skewer[i] <- sum((x[,i] - mx[i])^3,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^3)
            }
    }
    return(skewer)
}


#modified November 24, 2010 to use an unbiased estimator of kurtosis as the default
"kurtosi" <- 
function (x, na.rm = TRUE,type=1) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
       if (is.matrix(x) ) { mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
       
         sdx <- sd(x,na.rm=na.rm)
        kurt <- sum((x - mx)^4)/(length(x) *sd(x)^4)  -3
        } else {
    
    kurt <- rep(NA,dim(x)[2])
  #  mx <- mean(x,na.rm=na.rm)
  if (is.matrix(x) ) { mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
    
    sdx <- sd(x,na.rm=na.rm)
   
    for (i in 1:dim(x)[2]) {
    n <- length(x[!is.na(x[,i]),i])
   if(type !=1) { kurt[i] <- sum((x[,i] - mx[i])^4,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^4)  -3} else {
     xi <- x[,i]-mx[i]
     kurt[i] <- n*(n + 1)*sum((x[,i] - mx[i])^4,  na.rm = na.rm)/( (n - 1)*(n - 2)*(n - 3)*(sum((x[,i] - mx[i])^2,na.rm = na.rm)/(n - 1))^2)  -3 *(n- 1)^2 /((n - 2)*(n - 3)) } }
    names(kurt) <- colnames(x)
    }
    return(kurt)
}


#added November 15, 2010
#adapted from the mult.norm function of the QuantPsyc package

"mardia" <-
function(x,na.rm=TRUE,plot=TRUE) {
 cl <- match.call()
x <- as.matrix(x)     #in case it was a dataframe
if(na.rm) x <- na.omit(x)
n <- dim(x)[1]
p <- dim(x)[2]
x <- scale(x,scale=FALSE)  #zero center
S <- cov(x)
S.inv <- solve(S)
D <- x %*% S.inv %*% t(x)
b1p <- sum(D^3)/n^2
b2p <- tr(D^2)/n 
chi.df <- p*(p+1)*(p+2)/6
k <- (p+1)*(n+1)*(n+3)/(n*((n+1)*(p+1) -6))

small.skew <- n*k*b1p/6
M.skew <- n*b1p/6
M.kurt <- (b2p - p * (p+2))*sqrt(n/(8*p*(p+2)))
p.skew <- 1-pchisq(M.skew,chi.df)
p.small <- 1 - pchisq(small.skew,chi.df)
p.kurt <- 2*(1- pnorm(abs(M.kurt)))
d =sqrt(diag(D))
if(plot) {qqnorm(d)
          qqline(d)}
results <- list(n.obs=n,n.var=p, b1p = b1p,b2p = b2p,skew=M.skew,small.skew=small.skew,p.skew=p.skew,p.small=p.small,kurtosis=M.kurt,p.kurt=p.kurt,d = d,Call=cl)
class(results) <- c("psych","mardia")
return(results)
}


