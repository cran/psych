"scoreOverlap" <-
function(keys,r,correct=TRUE,SMC=TRUE,av.r=TRUE,item.smc=NULL,impute=TRUE) { #function to score clusters according to the key matrix, correcting for item overlap
				
 tol=sqrt(.Machine$double.eps)    #machine accuracy
 cl <- match.call()
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 if ((dim(r)[1] != dim(r)[2]) ) {r <- cor(r,use="pairwise")}
# if(any(is.na(r))) {SMC=FALSE
#                 warning("Missing values in the correlation matrix do not allow for SMC's to be found")}

 if(SMC && is.null(item.smc)) {item.smc <- smc(r)} else {item.smc <- rep(1,dim(r)[1]) }
                                   
 if(all(item.smc ==1)) SMC <- FALSE
 covar <- t(keys) %*% r %*% keys    #matrix algebra is our friend 
 var <- diag(covar)    #these are the scale variances
 n.keys <- ncol(keys)

 raw.r  <- cov2cor(covar)
 key.var <- diag(t(keys) %*% keys)
 key.smc <- t(keys) %*% item.smc  
 key.alpha <- ((var-key.var)/var)*(key.var/(key.var-1))
 key.lambda6 <-  (var - key.var + key.smc)/var
 key.alpha[is.nan(key.alpha)] <- 1           #if only 1 variable to the cluster, then alpha is undefined
 key.alpha[!is.finite(key.alpha)] <- 1   
 key.av.r <- key.alpha/(key.var - key.alpha*(key.var-1))  #alpha 1 = average r
 colnames(raw.r) <- rownames(raw.r)  <- colnames(keys)
 names(key.lambda6) <- colnames(keys)
 key.lambda6 <- drop(key.lambda6)
 
 n.keys <- ncol(keys)
 sn <- key.av.r * key.var/(1-key.av.r)
 item.cov <- t(keys) %*% r 
 
 raw.cov <- item.cov %*% keys
 adj.cov <- raw.cov
 
 #now adjust them
 
 for (i in 1:(n.keys)) {
    for (j in 1:i) {
 if(av.r) { adj.cov[i,j] <- adj.cov[j,i]<- raw.cov[i,j] - sum(keys[,i] * keys[,j] ) + sum(keys[,i] * keys[,j] *  sqrt(key.av.r[i] * key.av.r[j]))
  } else {
     adj.cov[i,j] <- adj.cov[j,i]<- raw.cov[i,j] - sum(keys[,i] * keys[,j] )+ sum( keys[,i] * keys[,j] * sqrt(item.smc[abs(keys[,i])>0]*item.smc[abs(keys[,j]) >0 ]))
 }
    } }

diag(adj.cov) <- diag(raw.cov)
adj.r <- cov2cor(adj.cov)


 
 
 
 if (correct) {cluster.corrected <- correct.cor(adj.r,t(key.alpha))
 result <- list(cor=adj.r,sd=sqrt(var),corrected= cluster.corrected,alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,Call=cl)
 }  #correct for attenuation
 else {
result <- list(cor=adj.r,sd=sqrt(var),alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,Call=cl)}
 class(result) <- c ("psych", "overlap")
 return(result)}
 #modified 01/11/15 to find r if not a square matrix
