"factor.fit" <-
function (r,f) {
     r2 <-sum( r*r)
     rstar <- factor.residuals(r,f)
     rstar2 <- sum(rstar*rstar)
     fit<- 1- rstar2/r2
     return(fit) }
     
