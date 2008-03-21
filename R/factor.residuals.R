"factor.residuals" <-
function(r, f) {
 if(is.matrix(f)) {
   rstar <- r - factor.model(f)} else {
   f <- f$loadings
    rstar <- r - factor.model(f)}
   return(rstar)}

