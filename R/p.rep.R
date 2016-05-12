"p.rep" <- function(p=.05,n=NULL,twotailed=FALSE) {
df <- n - 2
  if(twotailed) p <- 2*p
 p.rep  <- pnorm(qnorm((1-p))/sqrt(2))
 if (!is.null(n)) { t <- -qt(p/2,df)  
  r.equiv <-  sqrt(t^2/(t^2 + df))
  dprime = 2*t*sqrt(1/df)
  return(list(p.rep=p.rep,d.prime=dprime,r.equiv=r.equiv)) } else {
      
  return(p.rep)}
  }
  
"p.rep.f" <- function(F,df2,twotailed=FALSE) {
 p <- pf(F,1,df2,lower.tail=FALSE)

   dprime = sqrt(4*F/(df2))
 if(twotailed) p <- 2*p
 p.rep  <- pnorm(qnorm(1-p)/sqrt(2))
 if(twotailed) p <- 2*p
 r.equiv <- sqrt(F/(F+df2))
 return(list(p.rep=p.rep,dprime=dprime,prob=p,r.equiv=r.equiv))}
  
"p.rep.r" <- function(r,n,twotailed=TRUE) {
 dprime <- 2*r/sqrt(1-r^2)
 sigmad <- sqrt(4/(n-4))
 z <- dprime/sigmad
 p <- 1- pnorm(z)
 p.rep  <- pnorm(qnorm((1-p))/sqrt(2))
 if(twotailed) p <- 2*p

 return(list(p.rep=p.rep,dprime=dprime,prob=p))}
  
 "p.rep.t" <- function(t,df,df2=NULL,twotailed=TRUE) {
if (is.null(df2)) { 
     dprime = 2*t/sqrt(df)
 nc <- 1 } else { n1 <- df+1
           n2 <- df2+1
           df <- df + df2
           nc <- ((n1+n2)/2) / ((2*n1 * n2)/(n1+n2)) #average n /n harmonic
       dprime <- (2* t /sqrt(df)) * sqrt(nc)
       } 
 p <- pt(t,df,lower.tail=FALSE)
  r.equiv <-  sqrt(t^2/(t^2 + df))  
 if(twotailed) p <- 2*p
 p.rep  <- pnorm(qnorm((1-p))/sqrt(2))
  return(list(p.rep=p.rep,dprime=dprime,prob=p,r.equiv=r.equiv))}
  


"p.rep.z" <- function(z,n,twotailed=TRUE) {

     dprime = 2*z*sqrt(1/n)
 p <- (1-pnorm(z ))
  r.equiv <-  sqrt(z^2/(z^2 + n))
 if(twotailed) p <- 2*p
 p.rep  <- pnorm(qnorm((1-p))/sqrt(2))
  return(list(p.rep=p.rep,dprime=dprime,prob=p,r.equiv=r.equiv))}
  
