"scoreOverlap" <-
function(keys,r,correct=TRUE,SMC=TRUE,av.r=TRUE,item.smc=NULL,impute=TRUE,select=TRUE) { #function to score clusters according to the key matrix, correcting for item overlap
				
 tol=sqrt(.Machine$double.eps)    #machine accuracy
 cl <- match.call()
 bad <- FALSE
 key.list <- keys  #we need to keep this as a list for later  
 if(is.list(keys) & (!is.data.frame(keys))) { if (select) {
    select <- selectFromKeyslist(colnames(r),keys)
 # select <- sub("-","",unlist(keys))   #added April 7, 2017
      select <- select[!duplicated(select)]
      }  else {select <- 1:ncol(r) } 
      
 if (!isCorrelation(r,na.rm=TRUE)) {r <- cor(r[,select],use="pairwise")} else {r <- r[select,select]}
 keys <- make.keys(r,keys)}  #added 9/9/16    (and then modified March 4, 2017
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
 if ((dim(r)[1] != dim(r)[2]) ) {r <- cor(r,use="pairwise")}
 if(any(abs(r[!is.na(r)]) > 1)) warning("Something is seriously wrong with the correlation matrix, some correlations had absolute values > 1!  Please check your data.")
 if(any(is.na(r))) {
  #              SMC=FALSE
 #                warning("Missing values in the correlation matrix do not allow for SMC's to be found")
                 bad <- TRUE}

 if(SMC && is.null(item.smc)) {item.smc <- smc(r)} else {
         diag(r) <- NA
         item.smc <- apply(r,1,function(x) max(abs(x),na.rm=TRUE))
         item.smc[is.infinite(item.smc) ] <- 1 
         diag(r) <- 1}
                                   
 if(all(item.smc ==1)) SMC <- FALSE
 if(!bad) {covar <- t(keys) %*% r %*% keys} else  #matrix algebra is our friend 
     {#covar<- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))*x,na.rm=TRUE))  #matrix multiplication without matrices!
     covar <- score.na(keys,r,cor=FALSE)
  }

 var <- diag(covar)    #these are the scale variances
 n.keys <- ncol(keys)
 item.var <- item.smc
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
 
if(!bad) { item.cov <- t(keys) %*% r    #the normal case is to have all correlations
         raw.cov <- item.cov %*% keys} else {  
         item.cov <- apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))  #some correlations are NA have to adjust
         raw.cov <-  apply(keys,2,function(x) colSums(item.cov*x,na.rm=TRUE))
         item.cov <- t(item.cov)
   }
 adj.cov <- raw.cov 
 
 #now adjust them
 
  med.r <- rep(NA, n.keys)

 for (i in 1:(n.keys)) {
    temp <- keys[,i][abs(keys[,i]) > 0]
  temp <- diag(temp,nrow=length(temp))
   small.r <- r[abs(keys[,i])>0,abs(keys[,i])>0]
   #small.r <- temp %*% small.r %*% temp   #this is just flipping the signs, but will not work with missing data
   if(NROW(temp) > 1) small.r <- apply(temp,2, function(x) colSums(apply(temp,2, function(x) colSums(small.r * x,na.rm=TRUE))*x,na.rm=TRUE))


    med.r[i]  <- median(small.r[lower.tri(small.r)],na.rm=TRUE)  

    for (j in 1:i) {
   
 if(av.r) { adj.cov[i,j] <- adj.cov[j,i]<- raw.cov[i,j] - sum(keys[,i] * keys[,j] ) + sum(keys[,i] * keys[,j] *  sqrt(key.av.r[i] * key.av.r[j]))
  } else {
     adj.cov[i,j] <- adj.cov[j,i] <- raw.cov[i,j] - sum(keys[,i] * keys[,j] )+ sum( keys[,i] * keys[,j] * sqrt(item.smc[i]* abs(keys[,i])*item.smc[j]*abs(keys[,j]) ))
 
 }
    } }

scale.var <- diag(raw.cov)

diag(adj.cov) <- diag(raw.cov)
adj.r <- cov2cor(adj.cov)   #this is the overlap adjusted correlations

#find the MIMS values  (Average within cluster/scale items)
scale.size <- outer(key.var,key.var)
MIMS <- adj.cov/scale.size
diag(MIMS)<- key.av.r   



#adjust the item.cov for item overlap
#we do this by replacing the diagonal of the r matrix with the item.var (probably an smc, perhaps a maximum value)

diag(r) <- item.var
if(!bad) { item.cov <- t(keys) %*% r    #the normal case is to have all correlations
        } else {  
         item.cov <- t(apply(keys,2,function(x) colMeans(r*x,na.rm=TRUE)) *NROW(keys))  #some correlations are NA
         }

 if(n.keys > 1) {
    item.cor <-   sqrt(diag(1/(key.lambda6*scale.var))) %*% (item.cov)  # %*% diag(1/sqrt(item.var))
    rownames(item.cor) <- colnames(keys)
    } else {
      item.cor <- r %*% keys /sqrt(key.lambda6*scale.var) }
   colnames(item.cor) <- colnames(r)
   item.cor <- t(item.cor)
   names(med.r) <- colnames(keys)


#find the Multi-Item Multi Trait item x scale correlations
 MIMT <- matrix(NA,n.keys,n.keys)
 for (i in 1:(n.keys)) {
    temp <- keys[,i][abs(keys[,i]) > 0]
    flip.item <- temp * item.cor[names(temp),]
   if(length(names(temp)) > 1) { MIMT[i,] <- colMeans(flip.item[names(temp),])} else {MIMT[i,] <- flip.item}
    }
  colnames(MIMT) <- rownames(MIMT) <- colnames(keys)
 
 good <- scale_quality(adj.r,item.cor,key.list)
 names(good) <- names(key.list)
 
 if (correct) {cluster.corrected <- correct.cor(adj.r,t(key.alpha))
 result <- list(cor=adj.r,sd=sqrt(var),corrected= cluster.corrected,alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,item.cor=item.cor,med.r=med.r,quality=good,MIMS=MIMS,MIMT=MIMT,Call=cl)
 }  #correct for attenuation
 else {
result <- list(cor=adj.r,sd=sqrt(var),alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,item.cor=item.cor,med.r=med.r,Call=cl)}
 class(result) <- c ("psych", "overlap")
 return(result)}
 #modified 01/11/15 to find r if not a square matrix
#modifed 03/05/15 to do pseudo matrix multiplication in case of missing data 


scale_quality = function(phi,r,keys) { #switched from . to _ 6/20/23
nvar <- NROW(r)
nscale <- NCOL(r)
good <- rep(0,length(keys))
best <- apply(abs(r),1, which.max)
for(i in 1:length(keys)) {
select <- selectFromKeys(keys[i])
good [i] <-  sum(best[select] == i)
good[i] <- good[i]/length(select)
}

return(good)
}


# scale.quality = function(n.obs,phi,r,keys) {
# nvar <- NROW(r)
# nscale <- NCOL(r)
# best <- good <- rep(0,length(keys))
# best <- apply(abs(r),1, which.max)
# for(i in 1:length(keys)) {
# select <- selectFromKeys(keys[i])
# good [i] <-  sum(best[select] == i)
# }
# return(good)
# }

