"fisherz" <-
function(rho)  {0.5*log((1+rho)/(1-rho)) }   #converts r to z  

"fisherz2r" <-
  function(z) {(exp(2*z)-1)/(1+exp(2*z)) }   #converts back again
  
"r2d" <- 
function(rho) {2*rho/sqrt(1-rho^2)}

"d2r" <- 
function(d) {d/sqrt(d^2+4)}

"t2r" <- function(t,df) {t^2/(t^2 + df)}

"g2r" <- function(g,df,n) {g/sqrt(g^2 + 4*df/n)}

"chi2r" <- function(chi,n) {sqrt(chi/n)}

"r2chi" <- function(rho,n) { chi <- rho^2 *n}




  

  

