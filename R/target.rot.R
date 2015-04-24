"target.rot" <- 
function (x,keys=NULL) 
{
 if(!is.matrix(x) & !is.data.frame(x) )  {
        if(!is.null(x$loadings)) x <- as.matrix(x$loadings)
      } else {x <- x}   
    if (ncol(x) < 2) 
        return(x)
    dn <- dimnames(x)
    if(is.null(keys)) {Q <- factor2cluster(x)} else {Q <- keys}
    Q <- as.matrix(Q)
    if(dim(Q)[2] < 2) {stop("Cluster structure produces 1 cluster. Rotation is not meaningful with less than 2 factors")}
    U <- lm.fit(x, Q)$coefficients
    d <- diag(solve(t(U) %*% U))
    U <- U %*% diag(sqrt(d))
    dimnames(U) <- NULL
    z <- x %*% U
    
    ui <- solve(U)
    Phi <- ui %*% t(ui)
    dimnames(z) <- dn
    class(z) <- "loadings"
    result <- list(loadings = z, rotmat = U,Phi = Phi)
    class(result) <- c("psych","fa")
    return(result)
}
#Based upon Promax which was taken from the promax function with the addition of returning the angles between factors
#based upon a suggestion to the R-help news group by Ulrich Keller and John Fox. 
#if keys is null, this is the Promax function
#if keys are not null, this becomes a targeted rotation function similar to that suggested by Michael Brown
#created April 6, 2009 with the assistance of Pat Shrout and Steve Miller
#a better model is to call TargetQ


