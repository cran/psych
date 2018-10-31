"fisherz" <-
function(rho)  {0.5*log((1+rho)/(1-rho)) }   #converts r to z  

"fisherz2r" <-
  function(z) {(exp(2*z)-1)/(1+exp(2*z)) }   #converts back again
  
"r2d" <- 
function(rho) {2*rho/sqrt(1-rho^2)}

"d2r" <- 
function(d) {d/sqrt(d^2+4)}

#added sign correction October 8, 2018
"t2r" <- function(t,df) {sign(t) * sqrt(t^2/(t^2 + df))}  #fixed April 27, 2017

"g2r" <- function(g,df,n) {sign(g) * g/sqrt(g^2 + 4*df/n)}

"chi2r" <- function(chi2,n) {sqrt(chi2/n)}

"r2chi" <- function(rho,n) { chi2 <-( rho^2 *n)}

"cor2cov" <- "r2c" <- function(rho,sigma) { sigma <- diag(sigma)
cov <- sigma %*% rho %*% sigma
colnames(cov) <- rownames(cov) <- colnames(rho)
return(cov)}


  

  

