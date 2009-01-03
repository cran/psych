"ICC" <- 
  function(x,digits=2,alpha=.05) {
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
 
  #find the various F values from Shrout and Fleiss 
  F11 <- MSB/MSW
  df11n <- nc-1
  df11d <- nc*(nr-1)
  p11 <- 1-pf(F11,df11n,df11d)
  F21 <- MSB/MSE
  df21n <- nc-1
  df21d <- (nc-1)*(nr-1)
  p21 <- 1-pf(F21,df21n,df21d)
  F31 <- F21
  # results <- list(ICC1=round(ICC1,digits),ICC2=round(ICC2,digits),ICC3 = round(ICC3,digits),ICC12=round(ICC12,digits),ICC22=round(ICC22,digits),ICC32=round(ICC32,digits))
 # results <- t(results)
 
 results <- data.frame(matrix(NA,ncol=8,nrow=6))
 colnames(results ) <- c("type", "ICC","F","df1","df2","p","lower bound","upper bound")
 rownames(results) <- c("Single_raters_absolute","Single_random_raters","Single_fixed_raters", "Average_raters_absolute","Average_random_raters","Average_fixed_raters")
 results[1,1] = "ICC1"
 results[2,1] = "ICC2"
 results[3,1] = "ICC3"
 results[4,1] = "ICC1k"
 results[5,1] = "ICC2k"
 results[6,1] = "ICC3k"
 results[1,2] = ICC1
 results[2,2] = ICC2
 results[3,2] = ICC3
 results[4,2] = ICC12
 results[5,2] = ICC22
 results[6,2] = ICC32
 results[1,3] <- results[4,3] <- F11  
 results[2,3] <- F21
 results[3,3] <-  results[6,3] <-  results[5,3] <- F31 <- F21 
 results[5,3] <- F21  
 results[1,4] <-   results[4,4] <- df11n
 results[1,5] <-    results[4,5] <-df11d
 results[1,6] <- results[4,6] <- p11
 results[2,4] <-  results[3,4] <-  results[5,4] <-  results[6,4] <- df21n
 results[2,5] <-  results[3,5] <-  results[5,5] <-  results[6,5] <- df21d
 results[2,6]  <- results[5,6] <-  results[3,6]  <-results[6,6] <- p21
 
 #now find confidence limits
 #first, the easy ones
 F1L <- F11 / qf(1-alpha/2,df11n,df11d)  
 F1U <- F11 * qf(1-alpha/2,df11d,df11n)
 L1 <- (F1L-1)/(F1L+(nr-1))
 U1 <- (F1U -1)/(F1U+nr-1)
 F3L <- F31 /qf(1-alpha/2,df21n,df21d)
 F3U <- F31 * qf(1-alpha/2,df21d,df21n)
 results[1,7] <- L1
 results[1,8] <- U1
 results[3,7] <- (F3L-1)/(F3L+nr-1)
 results[3,8] <- (F3U-1)/(F3U+nr-1)
 results[4,7] <- 1- 1/F1L
 results[4,8] <- 1- 1/F1U
 results[6,7] <- 1- 1/F3L
 results[6,8] <- 1 - 1/F3U
 
 #the hard one is case 2   
 Fj <- MSJ/MSE
 vn <- (nr-1)*(nc-1)* ( (nr*ICC2*Fj+nc*(1+(nr-1)*ICC2) - nr*ICC2))^2
 vd <- (nc-1)*nr^2 * ICC2^2 * Fj^2 + (nc *(1 + (nr-1)*ICC2) - nr*ICC2)^2
 v <- vn/vd
 F3U <- qf(1-alpha/2,nc-1,v) 
 F3L <- qf(1-alpha/2,v,nc-1)
 
 L3 <- nc *(MSB- F3U*MSE)/(F3U*(nr*MSJ+(nr*nc-nr-nc)*MSE)+ nc*MSB)
 results[2,7] <- L3
 U3 <- nc *(F3L * MSB - MSE)/(nr * MSJ + (nr * nc - nr - nc)*MSE + nc * F3L * MSB)
  results[2,8] <- U3
 L3k <- L3 * nr/(1+ L3*(nr-1))
 U3k <- U3 * nr/(1+ U3*(nr-1))
 results[5,7] <- L3k
 results[5,8] <- U3k
 

 #clean up the output
 results[,2:8] <- round(results[,2:8],digits)
 
 
 return(results)
  
  }