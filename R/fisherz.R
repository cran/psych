"fisherz" <-
function(rho)  {0.5*log((1+rho)/(1-rho)) }   #converts r to z  

"fisherz2r" <-
  function(z) {(exp(2*z)-1)/(1+exp(2*z)) }   #converts back again
  
"r2d" <- 
function(rho) {2*rho/sqrt(1-rho^2)}

"d2r" <- 
function(d) {d/sqrt(d^2+4)}


  

  

