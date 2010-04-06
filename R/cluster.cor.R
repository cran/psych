"cluster.cor" <-
function(keys,r.mat,correct=TRUE,digits=2,SMC=TRUE) { #function to extract clusters according to the key vector
				#default is to correct for attenuation and show this above the diagonal
				#find the correlation matrix of scales made up of items defined in a keys matrix (e.g., extracted by factor2cluster) 
                #takes as input the keys matrix as well as a correlation matrix of all the items
 tol=sqrt(.Machine$double.eps)    #machine accuracy
 cl <- match.call()
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 r.mat[is.na(r.mat)] <- -9999999    #changes missing values to obviously incorrect values
 if(SMC) {item.smc <- smc(r.mat)} else {item.smc <- rep(1,dim(r.mat)[1])}
 covar <- t(keys) %*% r.mat %*% keys    #matrix algebra is our friend
 var <- diag(covar)    #these are the scale variances
 sd.inv <- 1/sqrt(var)
 ident.sd <- diag(sd.inv,ncol = length(sd.inv))
 cluster.correl <- ident.sd %*% covar  %*% ident.sd
 cluster.correl[abs(cluster.correl)  > (1+tol)] <- NA    #happens only if item correlations were missing -- we use 1+a little to avoid rounding problem
 key.var <- diag(t(keys) %*% keys)
 key.smc <- t(keys) %*% item.smc  
 key.alpha <- ((var-key.var)/var)*(key.var/(key.var-1))
 key.lambda6 <-  (var - key.var + key.smc)/var
 key.alpha[is.nan(key.alpha)] <- 1           #if only 1 variable to the cluster, then alpha is undefined
 key.alpha[!is.finite(key.alpha)] <- 1   
 key.av.r <- key.alpha/(key.var - key.alpha*(key.var-1))  #alpha 1 = average r
 colnames(cluster.correl) <- colnames(keys)
 rownames(cluster.correl) <- colnames(keys)
 names(key.lambda6) <- colnames(keys)
 key.lambda6 <- drop(key.lambda6)
 
# diag(r.mat) <- 0  
# row.range <- apply(r.mat,1,range,na.rm=TRUE)     
# row.max <- pmax(abs(row.range[1,]),abs(row.range[2,]))  #find the largest absolute similarity
 
 if (correct) {cluster.corrected <- correct.cor(cluster.correl,t(key.alpha))
 result <- list(cor=round(cluster.correl,digits),sd=round(sqrt(var),digits),corrected= round(cluster.corrected,digits),alpha=round(key.alpha,digits),av.r = round(key.av.r,digits),size=key.var,G6 =key.lambda6,Call=cl)
 }  #correct for attenuation
 else {
result <- list(cor=round(cluster.correl,digits),sd=round(sqrt(var),digits),alpha=key.alpha,av.r = round(key.av.r,2),size=key.var,G6 =key.lambda6,Call=cl)}
 class(result) <- c ("psych", "cluster.cor")
 return(result)}
#revised August 21, 2007 to add a smidgen to 1.0 in the looking for NAs.
#revised June 14, 2008 to add average.r
#revised August 25, 2009 to add lambda6
