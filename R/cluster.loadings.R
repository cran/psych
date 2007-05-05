cluster.loadings <- 
function (keys, r.mat, correct = TRUE,digits=2) 
{
    if (!is.matrix(keys)) {
        keys <- as.matrix(keys)}
      
    item.covar <- r.mat %*% keys          #item by cluster covariances
    covar <- t(keys) %*% item.covar  #variance/covariance of clusters
    var <- diag(covar)
    sd.inv <- 1/sqrt(var)
    
    key.count <- diag(t(keys) %*% keys)    #how many items in each cluster?
    if (correct) {
                  cluster.correct <- diag((key.count/(key.count - 1)))
                   for (i in 1:ncol(keys)) {
                         if (key.count[i]==1) {   #fix the case of 1 item keys
                   	     cluster.correct[i,i] <- 1
                   	} else { cluster.correct[i,i] <- key.count[i]/(key.count[i]-1)
                   	item.covar[,i] <- item.covar[,i] - keys[,i]} 
                   }   #i loop
                   correction.factor <- keys %*% cluster.correct
                   correction.factor[ correction.factor < 1] <- 1
                  item.covar <- item.covar * correction.factor
                  }
   
      
    ident.sd <- diag(sd.inv, ncol = length(sd.inv))
    c.loading <-  item.covar %*% ident.sd
    c.correl <- ident.sd %*% covar %*% ident.sd
    
    
   
    
    key.alpha <- ((var - key.count)/var) * (key.count/(key.count - 1))
    key.alpha[is.nan(key.alpha)] <- 1
    key.alpha[!is.finite(key.alpha)] <- 1
    colnames(c.loading) <- colnames(keys)
    colnames(c.correl) <- colnames(keys)
    rownames(c.correl) <- colnames(keys)
    rownames(c.loading) <- rownames(r.mat)
    
    if( ncol(keys) >1)  {cluster.corrected <- correct.cor(c.correl, t(key.alpha))} else {cluster.corrected <- c.correl}
    
    return(list(loadings=round(c.loading,digits), cor=round(c.correl,digits),corrected=round(cluster.corrected,digits), sd = round(sqrt(var),digits), alpha = round(key.alpha,digits),
            size = key.count))
    }
