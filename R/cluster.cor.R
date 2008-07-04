"cluster.cor" <-
function(keys,r.mat,correct=TRUE,digits=2) { #function to extract clusters according to the key vector
				#default is to correct for attenuation and show this above the diagonal
				#find the correlation matrix of scales made up of items defined in a keys matrix (e.g., extracted by factor2cluster) 
                #takes as input the keys matrix as well as a correlation matrix of all the items
 tol=sqrt(.Machine$double.eps)    #machine accuracy
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 r.mat[is.na(r.mat)] <- -9999999    #changes missing values to obviously incorrect values
 covar <- t(keys) %*% r.mat %*% keys    #matrix algebra is our friend
 var <- diag(covar)
 sd.inv <- 1/sqrt(var)
 ident.sd <- diag(sd.inv,ncol = length(sd.inv))
 cluster.correl <- ident.sd %*% covar  %*% ident.sd
 cluster.correl[abs(cluster.correl)  > (1+tol)] <- NA    #happens only if item correlations were missing -- we use 1+a little to avoid rounding problem
 key.var <- diag(t(keys) %*% keys)
 key.alpha <- ((var-key.var)/var)*(key.var/(key.var-1))
 key.alpha[is.nan(key.alpha)] <- 1           #if only 1 variable to the cluster, then alpha is undefined
 key.alpha[!is.finite(key.alpha)] <- 1   
 key.av.r <- key.alpha/(key.var - key.alpha*(key.var-1))  #alpha 1 = average r
 colnames(cluster.correl) <- colnames(keys)
 rownames(cluster.correl) <- colnames(keys)
 if (correct) {cluster.corrected <- correct.cor(cluster.correl,t(key.alpha))
 return(list(cor=round(cluster.correl,digits),sd=round(sqrt(var),digits),corrected= round(cluster.corrected,digits),alpha=round(key.alpha,digits),av.r = round(key.av.r,2),size=key.var))
 }  #correct for attenuation
 else {
 return(list(cor=round(cluster.correl,digits),sd=round(sqrt(var),digits),alpha=round(key.alpha,digits),av.r = round(key.av.r,2),size=key.var))}
 }
#revised August 21, 2007 to add a smidgen to 1.0 in the looking for NAs.
#revised June 14, 2008 to add average.r
