"tenberge" <- 
function(r,digits=2) {n <- dim(r)[2]
if(dim(r)[1] > n) {r <- cor(r,use="pairwise")}
vt <- sum(r)
off <- r
diag(off) <- 0
sum.off <- sum(off)
sumsq.off <- sum(off^2)
lambda1 <- n * sum(off)/((n-1)* vt)
lambda2 <- (sum.off+ sqrt(sumsq.off*n/(n-1)))/vt
lambda3 <- (sum.off +sqrt(sumsq.off+ sqrt((n * sum(off^4)/(n-1)))))/vt
lambda4  <- (sum.off +sqrt(sumsq.off+ sqrt(sum(off^4)+ sqrt((n * sum(off^8)/(n-1))))))/vt
return(list(mu0 = round(lambda1,digits),mu1=round(lambda2,digits),mu2 = round(lambda3,digits),mu3=round(lambda4,digits)))
}