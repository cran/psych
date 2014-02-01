"ICC" <- 
  function(x,missing=TRUE,alpha=.05) {
  cl <- match.call()
  if(is.matrix(x)) x <- data.frame(x)
  n.obs.original <- dim(x)[1]
  if(missing) { x1 <- try(na.fail(x))
             if(class(x1) == as.character("try-error"))  {
               x <- na.omit(x)
               n.obs <- dim(x)[1]
             stop("missing data were found for ",n.obs.original -n.obs, " cases \n Try again with na.omit  or set missing= FALSE and proceed at your own risk.")}}
  n.obs <- dim(x)[1]
  if(n.obs < n.obs.original) message("Warning, missing data were found for ",n.obs.original -n.obs, " cases")
  nj <- dim(x)[2]
  x.s <- stack(x)
  x.df <- data.frame(x.s,subs=rep(paste("S",1:n.obs,sep=""),nj))
  aov.x <- aov(values~subs+ind,data=x.df)
  s.aov <- summary(aov.x)
  stats <- matrix(unlist(s.aov),ncol=3,byrow=TRUE)
  MSB <- stats[3,1]
  MSW <- (stats[2,2] + stats[2,3])/(stats[1,2] + stats[1,3])
  MSJ <- stats[3,2]
  MSE <- stats[3,3]
  
  ICC1 <- (MSB- MSW)/(MSB+ (nj-1)*MSW)
  ICC2 <- (MSB- MSE)/(MSB + (nj-1)*MSE + nj*(MSJ-MSE)/n.obs)
  ICC3 <- (MSB - MSE)/(MSB+ (nj-1)*MSE)
  ICC12 <- (MSB-MSW)/(MSB)
  ICC22 <- (MSB- MSE)/(MSB +(MSJ-MSE)/n.obs)
  ICC32 <- (MSB-MSE)/MSB
 
  #find the various F values from Shrout and Fleiss 
  F11 <- MSB/MSW
  df11n <- n.obs-1
  df11d <- n.obs*(nj-1)
  p11 <- 1-pf(F11,df11n,df11d)
  F21 <- MSB/MSE
  df21n <- n.obs-1
  df21d <- (n.obs-1)*(nj-1)
  p21 <- 1-pf(F21,df21n,df21d)
  F31 <- F21
  
  
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
 #don't divide alpha level by 2  (changed on 2/1/14)
 F1L <- F11 / qf(1-alpha/2,df11n,df11d)  
 F1U <- F11 * qf(1-alpha/2,df11d,df11n)
 L1 <- (F1L-1)/(F1L+(nj-1))
 U1 <- (F1U -1)/(F1U+nj-1)
 F3L <- F31 /qf(1-alpha/2,df21n,df21d)
 F3U <- F31 * qf(1-alpha/2,df21d,df21n)
 results[1,7] <- L1
 results[1,8] <- U1
 results[3,7] <- (F3L-1)/(F3L+nj-1)
 results[3,8] <- (F3U-1)/(F3U+nj-1)
 results[4,7] <- 1- 1/F1L
 results[4,8] <- 1- 1/F1U
 results[6,7] <- 1- 1/F3L
 results[6,8] <- 1 - 1/F3U
 
 #the hard one is case 2   
 Fj <- MSJ/MSE
 vn <- (nj-1)*(n.obs-1)* ( (nj*ICC2*Fj+n.obs*(1+(nj-1)*ICC2) - nj*ICC2))^2
 vd <- (n.obs-1)*nj^2 * ICC2^2 * Fj^2 + (n.obs *(1 + (nj-1)*ICC2) - nj*ICC2)^2
 v <- vn/vd
 F3U <- qf(1-alpha/2,n.obs-1,v) 
 F3L <- qf(1-alpha/2,v,n.obs-1)
 
 L3 <- n.obs *(MSB- F3U*MSE)/(F3U*(nj*MSJ+(nj*n.obs-nj-n.obs)*MSE)+ n.obs*MSB)
 results[2,7] <- L3
 U3 <- n.obs *(F3L * MSB - MSE)/(nj * MSJ + (nj * n.obs - nj - n.obs)*MSE + n.obs * F3L * MSB)
  results[2,8] <- U3
 L3k <- L3 * nj/(1+ L3*(nj-1))
 U3k <- U3 * nj/(1+ U3*(nj-1))
 results[5,7] <- L3k
 results[5,8] <- U3k
 
 #clean up the output
 results[,2:8] <- results[,2:8]
 result <- list(results=results,summary=s.aov,stats=stats,MSW=MSW,Call=cl,n.obs=n.obs,n.judge=nj)
 class(result) <- c("psych","ICC")
 return(result)
  
  }