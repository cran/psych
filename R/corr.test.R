"corr.test" <-
function(x,y=NULL,use="pairwise",method="pearson",adjust="holm") {
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
p[p>1] <- 1
if (adjust !="none") {
  if (is.null(y)) {lp <- upper.tri(p)  #the case of a symmetric matrix
     pa <- p[lp]
     pa <- p.adjust(pa,adjust)
     p[upper.tri(p,diag=FALSE)] <- pa
  } else {
  p[] <- p.adjust(p,adjust)  #the case of an asymmetric matrix 
} }
result <- list(r = r,n=n,t=t,p=p,se=se,adjust=adjust,sym =sym, Call=cl)
class(result) <- c("psych", "corr.test")
return(result)
}
#modified 1/4/14 to report sample size once if they are all equal



"corr.p" <-
function(r,n,adjust="holm") {
cl <- match.call()
if(missing(n)) stop("The number of subjects must be specified")

t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
p[p>1] <- 1
if (adjust !="none") {
if(isSymmetric(p)) {
 lp <- upper.tri(p)  #the case of a symmetric matrix
     pa <- p[lp]
     pa <- p.adjust(pa,adjust)
     p[upper.tri(p,diag=FALSE)] <- pa
  }
  p[] <- p.adjust(p ,adjust)  #the case of an asymmetric matrix
} 
result <- list(r = r,n=n,t=t,p=p,Call=cl)
class(result) <- c("psych", "corr.test")
return(result)
}