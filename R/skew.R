#corrected May 7, 2007 
#modified October ,2011 to use apply for mean and sd
#modified April, 2012 to return 3 estimates, depending upon type
#partly based upon e1071  skewness and kurtosis
"skew" <- 
function (x, na.rm = TRUE,type=3) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
        sdx <- sd(x,na.rm=na.rm)
        mx <- mean(x)
        n <- length(x[!is.na(x)]) 
        switch(type,
        {skewer <- sqrt(n) *( sum((x - mx)^3,  na.rm = na.rm)/( sum((x - mx)^2,na.rm = na.rm)^(3/2)))}, #case 1
        {skewer <- n *sqrt(n-1) *( sum((x - mx)^3,  na.rm = na.rm)/((n-2) * sum((x - mx)^2,na.rm = na.rm)^(3/2)))}, #case 2
        {skewer <- sum((x - mx)^3)/(n * sd(x)^3) })  #case 3
        } else {
    
    skewer <- rep(NA,dim(x)[2])
    if (is.matrix(x)) {mx <- colMeans(x,na.rm=na.rm)} else {mx <- apply(x,2,mean,na.rm=na.rm)}
    sdx <- apply(x,2,sd,na.rm=na.rm)
    for (i in 1:dim(x)[2]) {
    n <- length(x[!is.na(x[,i]),i])
    switch(type,
    {skewer[i] <-sqrt(n) *( sum((x[,i] - mx[i])^3,  na.rm = na.rm)/( sum((x[,i] - mx[i])^2,na.rm = na.rm)^(3/2)))}, #type 1
    {skewer[i] <- n *sqrt(n-1) *( sum((x[,i] - mx[i])^3,  na.rm = na.rm)/((n-2) * sum((x[,i] - mx[i])^2,na.rm = na.rm)^(3/2)))},#type 2
    {skewer[i] <- sum((x[,i] - mx[i])^3,  na.rm = na.rm)/(n * sdx[i]^3)} #type 3
           ) #end switch
            } #end loop
    }
    return(skewer)
}


#modified November 24, 2010 to use an unbiased estimator of kurtosis as the default
#and again April 22, 2012 to include all three types 
"kurtosi" <- 
function (x, na.rm = TRUE,type=3) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
       if (is.matrix(x) ) { mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
        sdx <- sd(x,na.rm=na.rm)
        n <- length(x[!is.na(x)])
        switch(type,
        {kurt <- sum((x - mx)^4,  na.rm = na.rm)*n /(sum((x - mx)^2,na.rm = na.rm)^2)  -3},  #type 1
        { 
         kurt <- n*(n + 1)*sum((x - mx)^4,  na.rm = na.rm)/( (n - 1)*(n - 2)*(n - 3)*(sum((x - mx)^2,na.rm = na.rm)/(n - 1))^2)  -3 *(n- 1)^2 /((n - 2)*(n - 3)) }, # type 2
        {kurt <- sum((x - mx)^4)/(n *sdx^4)  -3} )  #	type 3
        } else {
    
    kurt <- rep(NA,dim(x)[2])
  #  mx <- mean(x,na.rm=na.rm)
   mx <-apply(x,2 ,mean,na.rm=na.rm)
  if(type==3)  sdx <- apply(x,2,sd,na.rm=na.rm)  
    
    for (i in 1:dim(x)[2]) {
    n <- length(x[!is.na(x[,i]),i])
    switch(type,
      { kurt[i] <- sum((x[,i] - mx[i])^4,  na.rm = na.rm)*length(x[,i]) /(sum((x[,i] - mx[i])^2,na.rm = na.rm)^2)  -3},  #type 1
     {
     xi <- x[,i]-mx[i]
     kurt[i] <- n*(n + 1)*sum((x[,i] - mx[i])^4,  na.rm = na.rm)/( (n - 1)*(n - 2)*(n - 3)*(sum((x[,i] - mx[i])^2,na.rm = na.rm)/(n - 1))^2)  -3 *(n- 1)^2 /((n - 2)*(n - 3)) }  #type 2
     ,
   {
    kurt[i] <- sum((x[,i] - mx[i])^4,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^4)  -3},  #type 3
    {NULL}) 
    names(kurt) <- colnames(x)
    }}
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


