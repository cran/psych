cluster.loadings <- 
function (keys, r.mat, correct = TRUE,SMC=TRUE) 
{  cl <- match.call()
    if (!is.matrix(keys)) {
        keys <- as.matrix(keys)}
    r.mat[is.na(r.mat)] <- -9999999 
    
 
   
    item.sd <- sqrt(diag(r.mat))
    item.covar <- r.mat %*% keys          #item by cluster covariances
    covar <- t(keys) %*% item.covar  #variance/covariance of clusters
   
    var <- diag(covar)
    sd.inv <- 1/sqrt(var)
    
    
  #items corrected for communality lead to the Guttman G6 estimate
    if(SMC)  {r.smc <- smc(r.mat)
              r.smc[r.smc < 0 ] <- 1  #for a very weird condition
            diag(r.mat) <- r.smc } else {
            diag(r.mat) <- 0
            item.max <- apply(r.mat,1,max)
            diag(r.mat) <- item.max}
    c.item.covar <- r.mat %*% keys #this uses the communality estimate and thus corrects for item overlap
    c.covar <- t(keys) %*% c.item.covar 
    c.var <- diag(c.covar)
    G6 <- c.var/var
    n.keys <- dim(keys)[2] 
    if(n.keys >1) {
    c.item.cor <- c.item.covar %*% sqrt(diag(1/c.var))/item.sd } else {c.item.cor <- c.item.covar/sqrt(c.var*item.sd) }
    
  
     
    key.count <- diag(t(keys) %*% keys)    #how many items in each cluster?
    if (correct) {
                  cluster.correct <- diag((key.count/(key.count - 1)))
                   for (i in 1:dim(keys)[2]) {
                         if (key.count[i]<2 ) {   #fix the case of 1 item keys
                   	     cluster.correct[i,i] <- 1
                   	} else { cluster.correct[i,i] <- key.count[i]/(key.count[i]-1)
                   	item.covar[,i] <- item.covar[,i] - keys[,i]} #subtract the variance of the item
                   }   #i loop
                   correction.factor <- keys %*% cluster.correct #put back average correlation for the item if it loads on the key
                   correction.factor[ correction.factor < 1] <- 1
                  item.covar <- item.covar * correction.factor
                  }
   
      
    ident.sd <- diag(sd.inv, ncol = length(sd.inv))
   
    c.correl <- ident.sd %*% covar %*% ident.sd
    p.loading <- try(c.item.cor %*% solve(c.correl))
    if(class(p.loading)=="try-error") {message('the correlation matrix was singular, pattern loadings not found, proceed with caution')
                         p.loading  <- c.item.cor}
    
     c.item.cor[abs(c.item.cor)  > 99999] <- NA
     c.correl[abs(c.correl) > 99999] <- NA
    
   
    
    key.alpha <- ((var - key.count)/var) * (key.count/(key.count - 1))
    key.alpha[is.nan(key.alpha)] <- 1
    key.alpha[!is.finite(key.alpha)] <- 1
    key.av.r <- key.alpha/(key.count - key.alpha*(key.count-1))  #alpha 1 = average r
    colnames(c.item.cor) <- colnames(keys)
    colnames(p.loading) <- colnames(keys)
    colnames(c.correl) <- colnames(keys)
    rownames(c.correl) <- colnames(keys)
    rownames(c.item.cor) <- rownames(r.mat)
    
    if( ncol(keys) >1)  {cluster.corrected <- correct.cor(c.correl, t(key.alpha))} else {cluster.corrected <- c.correl}
    
    results <- list(loadings=c.item.cor,pattern=p.loading, cor=c.correl,corrected=cluster.corrected, sd = sqrt(var), alpha = key.alpha,av.r = key.av.r,
            size = key.count,G6=G6,Call=cl)
    class(results) <- c("psych","cluster.loadings")
    return(results)
    }
