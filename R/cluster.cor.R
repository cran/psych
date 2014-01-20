"cluster.cor" <-

#added  smc.items   30.12/13 to reduce the number of calls to SMC to once
function(keys,r.mat,correct=TRUE,SMC=TRUE,item.smc=NULL,impute=TRUE) { #function to extract clusters according to the key vector
				#default is to correct for attenuation and show this above the diagonal
				#find the correlation matrix of scales made up of items defined in a keys matrix (e.g., extracted by factor2cluster) 
                #takes as input the keys matrix as well as a correlation matrix of all the items
 tol=sqrt(.Machine$double.eps)    #machine accuracy
 cl <- match.call()
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 if(any(is.na(r.mat))) {SMC=FALSE
                 warning("Missing values in the correlation matrix do not allow for SMC's to be found")}
 r.mat[is.na(r.mat)] <- -9999999    #changes missing values to obviously incorrect values
 if(SMC && is.null(item.smc)) {item.smc <- smc(r.mat)} else {item.smc <- rep(1,dim(r.mat)[1])}   #now done once in the main function
 
 covar <- t(keys) %*% r.mat %*% keys    #matrix algebra is our friend but slow if we are doing this in iclust
   #probably better to just modify the matrix for those two rows and columns that are changing if doing an iclust
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
#now, try to figure out the imputed correlation for the case of NAs.
if(any(is.na(cluster.correl)) && impute) { #find the missing values based upon average covariances rather than totals
#first, change those bad values back to NA
  warning('Some of the correlations were NA and were imputed')
  r.mat[r.mat < -1] <- NA
  n.keys <- ncol(keys)
  keys[keys==0] <- NA  #this will allow us to find average responses
  for (i in 1:n.keys) { 
     #first find which variables are screwed up
     if(any(is.na(cluster.correl[i,]))) {#fix them
        for (j in 1:n.keys) {if(is.na(cluster.correl[i,j])) {#fix it
        temp <- mean(colMeans((keys[,i] * r.mat),na.rm=TRUE) * keys[,j],na.rm=TRUE)*key.var[i]*key.var[j]  #this is the  average covariance times the number of items scored
        adjusted.r <- temp * ident.sd[i,i]* ident.sd[j,j]
        cluster.correl[i,j] <- adjusted.r
        }     
         }
  }
  } 
 }
 sn <- key.av.r * key.var/(1-key.av.r)
 if (correct) {cluster.corrected <- correct.cor(cluster.correl,t(key.alpha))
 result <- list(cor=cluster.correl,sd=sqrt(var),corrected= cluster.corrected,alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,Call=cl)
 }  #correct for attenuation
 else {
result <- list(cor=cluster.correl,sd=sqrt(var),alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,Call=cl)}
 class(result) <- c ("psych", "cluster.cor")
 return(result)}
#revised August 21, 2007 to add a smidgen to 1.0 in the looking for NAs.
#revised June 14, 2008 to add average.r
#revised August 25, 2009 to add lambda6
#revised December 2011 to remove digits -- this is all handled in the print function
#revised January 2012 to estimate values when we have missing values in the correlations
#74% of the time is spent doing matrix multiplication  70 seconds for a 400 x 400 problem.


"icluster.cor" <-
#added to speed up iclust by just combining rows

function(keys,r.mat,ivar,jvar,correct=TRUE,SMC=TRUE,item.smc=NULL,impute=TRUE) { #function to extract clusters according to the key vector
				#default is to correct for attenuation and show this above the diagonal
				#find the correlation matrix of scales made up of items defined in a keys matrix (e.g., extracted by factor2cluster) 
                #takes as input the keys matrix as well as a correlation matrix of all the items
 tol=sqrt(.Machine$double.eps)    #machine accuracy
 cl <- match.call()
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 
 #covar <- t(keys) %*% r.mat %*% keys    #matrix algebra is our friend but slow if we are doing this in iclust
covar <- r.mat[-jvar,-jvar]   #drop the jth column and row
covar[ivar,] <- covar[,ivar] <- r.mat[ivar,] + r.mat[jvar,]
   #probably better to just modify the matrix for those two rows and columns that are changing if doing an iclust
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
#now, try to figure out the imputed correlation for the case of NAs.
if(any(is.na(cluster.correl)) && impute) { #find the missing values based upon average covariances rather than totals
#first, change those bad values back to NA
  warning('Some of the correlations were NA and were imputed')
  r.mat[r.mat < -1] <- NA
  n.keys <- ncol(keys)
  keys[keys==0] <- NA  #this will allow us to find average responses
  for (i in 1:n.keys) { 
     #first find which variables are screwed up
     if(any(is.na(cluster.correl[i,]))) {#fix them
        for (j in 1:n.keys) {if(is.na(cluster.correl[i,j])) {#fix it
        temp <- mean(colMeans((keys[,i] * r.mat),na.rm=TRUE) * keys[,j],na.rm=TRUE)*key.var[i]*key.var[j]  #this is the  average covariance times the number of items scored
        adjusted.r <- temp * ident.sd[i,i]* ident.sd[j,j]
        cluster.correl[i,j] <- adjusted.r
        }     
         }
  }
  } 
 }
 if (correct) {cluster.corrected <- correct.cor(cluster.correl,t(key.alpha))
 result <- list(cor=cluster.correl,sd=sqrt(var),corrected= cluster.corrected,alpha=key.alpha,av.r = key.av.r,size=key.var,G6 =key.lambda6,Call=cl)
 }  #correct for attenuation
 else {
result <- list(cor=cluster.correl,sd=sqrt(var),alpha=key.alpha,av.r = key.av.r,size=key.var,G6 =key.lambda6,Call=cl)}
 class(result) <- c ("psych", "cluster.cor")
 return(result)}
#revised August 21, 2007 to add a smidgen to 1.0 in the looking for NAs.
#revised June 14, 2008 to add average.r
#revised August 25, 2009 to add lambda6
#revised December 2011 to remove digits -- this is all handled in the print function
#revised January 2012 to estimate values when we have missing values in the correlations