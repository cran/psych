"circ.tests" <-
function(loads,loading=TRUE,sorting=TRUE) {
 cl <- match.call()

circ.gap <- function(loads,loading=TRUE,sorting=TRUE) {
if (loading) {l <- loads$loadings} else { 
       l <- loads}
       l<- l[,1:2]
    commun=rowSums(l*l) 
    theta=sign(l[,2])*acos(l[,1]/sqrt(commun)) #vector angle in radians
    if(sorting) {theta<- sort(theta)}
    gaps <- diff(theta)
    test <- var(gaps)
    return(test)
    }


circ.fisher <- function(loads,loading=TRUE) {
if (loading) {l <- loads$loadings} else { 
       l <- loads}
       l<- l[,1:2]
   radius <- rowSums(l * l)
   test <- sd(radius)/mean(radius)
  return (test)
  }
  
  
  circ.rt <- function(loads,loading=TRUE) {
 if (loading) {l <- loads$loadings} else { 
       l <- loads}
       l<- l[,1:2]
       qmc <- rep(0,10)
       for (i in 0:9) {theta <- 5*i
       	rl <- factor.rotate(l,theta,1,2)
       	 l2 <- rl*rl
       qmc[i] <- sum(apply(l2,1,var)) }
       test <- sd(qmc)/mean(qmc)
 }
 
 circ.v2 <- function(loads,loading=TRUE) {
if (loading) {l <- loads$loadings} else { 
       l <- loads}
       l<- l[,1:2]
   crit <- rep(0,10)
       for (i in 0:9) {
       		theta <- 5*i
       		rl <- factor.rotate(l,theta,1,2)
       	 	l2 <- rl*rl
       	 	suml2 <- sum(l2)
       		crit[i] <- var(l2[,1]/suml2)
       }
       test <- sd(crit)/mean(crit)
  return (test)
  }
 

   gap.test <- circ.gap(loads,loading,sorting)
   fisher.test <- circ.fisher(loads,loading)
   rotation.test <- circ.rt(loads,loading)
   variance.test <- circ.v2(loads,loading)
   circ.tests <- list(gaps=gap.test,fisher=fisher.test,RT=rotation.test,VT=variance.test,Call=cl)
   class(circ.tests) <- c("psych","circ")
   return(circ.tests)
}