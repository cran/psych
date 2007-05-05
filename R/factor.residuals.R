"factor.residuals" <-
function(r, f) {
   rstar<- r- factor.model(f)
   return(rstar)}

