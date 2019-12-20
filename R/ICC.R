#revised March, 9, 2018 to include lmer option
#this makes it much faster for large data sets, particularly with missing data

"ICC" <- 
  function(x,missing=TRUE,alpha=.05,lmer=TRUE,check.keys=FALSE) {
  cl <- match.call()
  if(is.matrix(x)) x <- data.frame(x)
  n.obs.original <- dim(x)[1]
 
  if(missing &!lmer) { x1 <- try(na.fail(x))
             if(inherits(x1, as.character("try-error")))  {
               x <- na.omit(x)
               n.obs <- dim(x)[1]
             stop("missing data were found for ",n.obs.original -n.obs, " cases \n Try again with na.omit  or set missing= FALSE and proceed at your own risk, or try the lmer=TRUE option")}}
  n.obs <- dim(x)[1]
  if(n.obs < n.obs.original) message("Warning, missing data were found for ",n.obs.original -n.obs, " cases")
  # n.obs <- NROW(x)
 items <- colnames(x)
n.items  <- NCOL(x)

if(check.keys) {
  min.item <- min(x[items],na.rm=TRUE)
  max.item <- max(x[items],na.rm=TRUE)
  
  p1 <- pca(x)
  keys <- rep(1,n.items)  
  if(any(p1$loadings < 0)) {
   message("Some items were negatively correlated with total scale and were automatically reversed.\n This is indicated by a negative sign for the variable name.") 
     keys[p1$loadings < 0] <- -1 }
 if(any(keys < 0)) {   
# #then find the x and y scores
 newx <- t(t(x) * keys + (keys < 0)*  (max.item + min.item) )  #these are now rescaled in the keyed direction 
 names(keys) <-colnames(x)
 cat("\ reversed items\n  ",names(keys)[keys < 0])
 x <- as.data.frame(newx)
 }  else {  message("Some items were negatively correlated with total scale.  To correct this problem, run again with check.keys=TRUE") }
 }
# 
  
  
  nj <- dim(x)[2]
  x.s <- stack(x)
  x.df <- data.frame(x.s,subs=rep(paste("S",1:n.obs,sep=""),nj))
  #choose to either do lmer or aov
   if(lmer) {
#lmer
colnames(x.df ) <- c("values","items","id")    #this makes it simpler to understand
    mod.lmer <- lme4::lmer(values ~ 1 + (1 | id) + (1 | items),
  data=x.df,na.action=na.omit)
  vc <- lme4::VarCorr(mod.lmer)
  
  MS_id <- vc$id[1,1]
  MS_items <- vc$items[1,1]
  MSE <- error <- MS_resid <- (attributes(vc)$sc)^2 
  MS.df <- data.frame(variance= c(MS_id ,MS_items, MS_resid,NA))
  rownames(MS.df) <- c("ID","Items", "Residual","Total") 

 MS.df["Total",]  <- sum(MS.df[1:3,1],na.rm=TRUE)
 MS.df["Percent"] <- MS.df/MS.df["Total",1]
lmer.MS <- MS.df  #save these
#convert to AOV equivalents
 MSB <- nj * MS_id + error
 MSJ <- n.obs * MS_items + error
 MSW <- error + MS_items
 stats <- matrix(NA,ncol=3,nrow=5)  #create an anova equivalent table
stats[1,1] <- dfB <-  n.obs -1
stats[1,2] <- dfJ <-  nj - 1
stats[1,3] <- dfE <-  ( n.obs -1) * ( nj - 1)
stats[2,1] <- MSB *( n.obs -1)
stats[2,2] <- MSJ *(nj -1)
stats[2,3] <- MSE * (n.obs-1) *(nj-1)
 stats[3,1] <- MSB
 stats[3,2] <- MSJ
 stats[3,3] <- MSE
 stats[4,1] <- FB <-  MSB/MSE
 stats[4,2] <- FJ <-  MSJ/MSE
 stats[5,1] <-  -expm1(pf(FB,dfB,dfE,log.p=TRUE))  
 stats[5,2] <- -expm1(pf(FJ,dfJ,dfE,log.p=TRUE))  
# s.aov  <-  MS.df
s.aov <- mod.lmer
 
} else {
#AOV
  aov.x <- aov(values~subs+ind,data=x.df)
  s.aov <- summary(aov.x)
  stats <- matrix(unlist(s.aov),ncol=3,byrow=TRUE)
  MSB <- stats[3,1]
  MSW <- (stats[2,2] + stats[2,3])/(stats[1,2] + stats[1,3])
  MSJ <- stats[3,2]
  MSE <- stats[3,3]
 MS.df <- NULL
 }
 colnames(stats) <- c("subjects","Judges", "Residual")
 rownames(stats) <- c("df","SumSq","MS","F","p")
  
  
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
 # p11 <- 1-pf(F11,df11n,df11d)
   p11 <-  -expm1(pf(F11,df11n,df11d,log.p=TRUE))  
  F21 <- MSB/MSE
  df21n <- n.obs-1
  df21d <- (n.obs-1)*(nj-1)
 # p21 <- 1-pf(F21,df21n,df21d)
   p21 <-  - expm1(pf(F21,df21n,df21d,log.p=TRUE)) 
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
 #fixed again? onf 5/21/19
  F1L <- F11 / qf(1-alpha,df11n,df11d)  
 F1U <- F11 * qf(1-alpha,df11d,df11n)
 L1 <- (F1L-1)/(F1L+(nj-1))
 U1 <- (F1U -1)/(F1U+nj-1)
 F3L <- F31 /qf(1-alpha,df21n,df21d)
 F3U <- F31 * qf(1-alpha,df21d,df21n)
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
 F3U <- qf(1-alpha,n.obs-1,v) 
 F3L <- qf(1-alpha,v,n.obs-1)
 
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
 result <- list(results=results,summary=s.aov,stats=stats,MSW=MSW,lme = MS.df,Call=cl,n.obs=n.obs,n.judge=nj)
 class(result) <- c("psych","ICC")
 return(result)
  }
  
  
