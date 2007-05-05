"eigen.loadings" <-
function (x) { #convert eigen vectors to loadings by unnormalizing them
               #used if using princomp or princ, not needed for principal
    n <- length(x$values)
    x$values[ x$values<0]<- 0
    fix<-sqrt(x$values)
    result<- x$vectors * rep(fix, each = n)}

