"fa.extension" <-
  function(Ro,Roe,fo,correct=TRUE) {
 cl <- match.call()
  w <- fo$weights
  fo <- fo$loadings
  r2 <- t(w) %*% fo
  fe <- t(t(w) %*% Roe)
  d <-diag(t(fo) %*% w)
 if(correct) fe <- (fe/d)
 colnames(fe) <- colnames(fo)
 result <- list(loadings = fe,Call=cl)
 class(result) <- c("psych","extension")
 return(result)
}