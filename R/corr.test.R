"corr.test" <-
function(x,y=NULL,use="pairwise",method="pearson",adjust="holm",alpha=.05){
cl <- match.call()
if(is.null(y)) {r <- cor(x,use=use,method=method)
 sym <- TRUE
n <- t(!is.na(x)) %*% (!is.na(x))
} else {r <- cor(x,y,use=use,method=method)
  sym=FALSE
n <- t(!is.na(x)) %*% (!is.na(y))}
if((use=="complete") | (min(n) == max(n))) n <- min(n)
t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
se <- sqrt((1-r*r)/(n-2))

nvar <- ncol(r)

p[p>1] <- 1
if (adjust !="none") {
  if (is.null(y)) {lp <- upper.tri(p)  #the case of a symmetric matrix
     pa <- p[lp]
     pa <- p.adjust(pa,adjust)
     p[upper.tri(p,diag=FALSE)] <- pa
  } else {
  p[] <- p.adjust(p,adjust)  #the case of an asymmetric matrix 
} }

  z <- fisherz(r[lower.tri(r)])
  
   if (min(n) < 4) {
      warning("Number of subjects must be greater than 3 to find confidence intervals.")
   }
    
    se <- 1/sqrt(n[lower.tri(n)] - 3)
    alpha <- 1-alpha/2
    dif <- qnorm(alpha)
    lower <- fisherz2r(z - dif * se)
    upper <- fisherz2r(z + dif * se)
  ci <- data.frame(lower=lower,r=r[lower.tri(r)],upper=upper,p=p[lower.tri(p)])
      cnR <- abbreviate(colnames(r),minlength=5) 
      k <- 1
      for(i in 2:nvar) {for (j in 1:(i-1)) {
      rownames(ci)[k] <- paste(cnR[j],cnR[i],sep="-")
      k<- k +1 }}

result <- list(r = r,n=n,t=t,p=p,se=se,adjust=adjust,sym =sym,ci=ci, Call=cl)
class(result) <- c("psych", "corr.test")
return(result)
}
#modified 1/4/14 to report sample size once if they are all equal
#modified 3/12/14 to report conficence intervals (suggested by Alexander Weiss)



"corr.p" <-
function(r,n,adjust="holm",alpha=.05) {
cl <- match.call()
if(missing(n)) stop("The number of subjects must be specified")
sym <- FALSE
t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
p[p>1] <- 1
if (adjust !="none") {
if(isSymmetric(unclass(p))) {sym <- TRUE
 lp <- upper.tri(p)  #the case of a symmetric matrix
     pa <- p[lp]
     pa <- p.adjust(pa,adjust)
     p[upper.tri(p,diag=FALSE)] <- pa
  } else {
  p[] <- p.adjust(p ,adjust)  #the case of an asymmetric matrix
  sym <- FALSE}
} 
nvar <- ncol(r)
z <- fisherz(r[lower.tri(r)])
 if (min(n) < 4) {
      warning("Number of subjects must be greater than 3 to find confidence intervals.")
   }
   if(is.matrix(n)) {
   se <- 1/sqrt(n[lower.tri(n)] - 3) } else { se <- 1/sqrt(n - 3)}
    alpha <- 1-alpha/2
    dif <- qnorm(alpha)
    lower <- fisherz2r(z - dif * se)
    upper <- fisherz2r(z + dif * se)
  ci <- data.frame(lower=lower,r=r[lower.tri(r)],upper=upper,p=p[lower.tri(p)])
      cnR <- abbreviate(colnames(r),minlength=5) 
      k <- 1
      for(i in 2:nvar) {for (j in 1:(i-1)) {
      rownames(ci)[k] <- paste(cnR[j],cnR[i],sep="-")
      k<- k +1 }}
result <- list(r = r,n=n,t=t,p=p,sym=sym,adjust=adjust,ci=ci,Call=cl)
class(result) <- c("psych", "corr.p")
return(result)
}