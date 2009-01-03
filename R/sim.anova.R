sim.anova <- function(n=12,es1=1,es2=1,es3=1) {
IV1 <- c(rep(-1,n/2),rep(1,n/2))
IV2 <- rep(c(-1,1),n/2)
 
 y <- es1*IV1 + es2*IV2 + es3*IV1*IV2 + rnorm(n)
 y.df <- data.frame(y,IV1,IV2)
 }