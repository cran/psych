"factor.residuals" <-
function(r, f) {
 if(is.matrix(f)) {
   rstar <- r - factor.model(f)} else {
   Phi <- f$Phi
   f <- f$loadings
    rstar <- r - factor.model(f,Phi=Phi)}
   return(rstar)}


