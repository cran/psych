"ICC" <- 
  function(x,digits=2) {
  if(is.matrix(x)) x <- data.frame(x)
  nc <- dim(x)[1]
  nr <- dim(x)[2]
  x.s <- stack(x)
  x.df <- data.frame(x.s,subs=rep(paste("S",1:nc,sep=""),nr))
  aov.x <- aov(values~subs+ind,data=x.df)
  s.aov <- summary(aov.x)
  stats <- matrix(unlist(s.aov),ncol=3,byrow=TRUE)
  MSB <- stats[3,1]
  MSW <- (stats[2,2] + stats[2,3])/(stats[1,2] + stats[1,3])
  MSJ <- stats[3,2]
  MSE <- stats[3,3]
  
  ICC1 <- (MSB- MSW)/(MSB+ (nr-1)*MSW)
  ICC2 <- (MSB- MSE)/(MSB + (nr-1)*MSE + nr*(MSJ-MSE)/nc)
  ICC3 <- (MSB - MSE)/(MSB+ (nr-1)*MSE)
  ICC12 <- (MSB-MSW)/(MSB)
  ICC22 <- (MSB- MSE)/(MSB +(MSJ-MSE)/nc)
  ICC32 <- (MSB-MSE)/MSB
  results <- list(ICC1=round(ICC1,digits),ICC2=round(ICC2,digits),ICC3 = round(ICC3,digits),ICC12=round(ICC12,digits),ICC22=round(ICC22,digits),ICC32=round(ICC32,digits))
  results <- t(results)
  return(results)
  
  }