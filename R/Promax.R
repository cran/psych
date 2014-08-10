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
#based upon a suggestion to the R-help news group by Ulrich Keller and John Fox. 

#added May 31st following suggestions to R-Help by Gunter Nickel
"equamax" <- function(L, Tmat=diag(ncol(L)),  eps=1e-5, maxit=1000) {
kappa=ncol(L)/(2*nrow(L))
cfT(L, Tmat=diag(ncol(L)),  eps=eps, maxit=maxit)}


#based completely on the GPArotation  GPForth function
#modified to call the varimin function which is derived from the varimax function

varimin <- function(L, Tmat = diag(ncol(L)), normalize = FALSE, eps = 1e-05, 
    maxit = 1000) {
    GPForth(A=L,Tmat = diag(ncol(L)), normalize = normalize, eps = eps, 
    maxit = maxit, method = "varimin") }
    
 vgQ.varimin <- 
function (L) 
{
    QL <- sweep(L^2, 2, colMeans(L^2), "-")
    list(Gq = L * QL, f = sqrt(sum(diag(crossprod(QL))))^2/4, 
        Method = "varimin")
}
   

