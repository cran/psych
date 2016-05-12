"r.con" <- 
  function(rho,n,p=.95,twotailed=TRUE) {
   z <- fisherz(rho)
   if(n<4) {stop("number of subjects must be greater than 3")}
   se <- 1/sqrt(n-3)
   p <- 1-p 
   if(twotailed) p<- p/2
   dif <- qnorm(p)
   zlow <- z + dif*se
   zhigh <- z - dif*se
   ci <- c(zlow,zhigh)
   ci <- fisherz2r(ci)
   return(ci)
   }
 
 
   
"r2t" <- 
   function(rho,n) {
   return( rho*sqrt((n-2)/(1-rho^2))) }
   