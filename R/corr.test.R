"corr.test" <-
function(x,y=NULL,use="pairwise",method="pearson") {
cl <- match.call()
if(is.null(y)) {r <- cor(x,use=use,method=method)
n <- t(!is.na(x)) %*% (!is.na(x))
} else {r <- cor(x,y,use=use,method=method)
n <- t(!is.na(x)) %*% (!is.na(y))}
if(use=="complete") n <- min(n)
t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
p[p>1] <- 1

result <- list(r = r,n=n,t=t,p=p,Call=cl)
class(result) <- c("psych", "corr.test")
return(result)
}