"Promax" <- 
function (x, m = 4) 
{
 if(!is.matrix(x) & !is.data.frame(x) )  {
        if(!is.null(x$loadings)) x <- as.matrix(x$loadings)
      } else {x <- x}   
    if (ncol(x) < 2) 
        return(x)
    dn <- dimnames(x)
    xx <- varimax(x)
    x <- xx$loadings
    Q <- x * abs(x)^(m - 1)
    U <- lm.fit(x, Q)$coefficients
    d <- try(diag(solve(t(U) %*% U)),silent=TRUE)
     if(class(d)=="try-error") {warning("Factors are exactly uncorrelated and the model produces a singular matrix. An approximation is used")
        ev <- eigen(t(U) %*% U)
        ev$values[ev$values < .Machine$double.eps] <- 100 * .Machine$double.eps
        UU <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
       diag(UU)  <- 1
       d <- diag(solve(UU))}
    U <- U %*% diag(sqrt(d))
    dimnames(U) <- NULL
    z <- x %*% U
    U <- xx$rotmat %*% U
    ui <- solve(U)
    Phi <- ui %*% t(ui)
    dimnames(z) <- dn
    class(z) <- "loadings"
    result <- list(loadings = z, rotmat = U,Phi = Phi)
    class(result) <- c("psych","fa")
    return(result)
}
#obviously a direct copy of the promax function, with the addition of returning the angles between factors
#based upon a suggestion to the R-help news group by ulrich keller and John Fox. 