#Find scores from correlations for each subject
#this is a version scoreOverlap adjusted the case of many subjects
#impute missing values?
#The correlations for each subject are found by statsBy
#and are stored as a list of matrices in the r element
#The trick is how to handle missing correlations 
#Various fixes 1/15/22
"scoreBy" <-
function(keys,stats,correct=TRUE,SMC=TRUE,av.r=TRUE,item.smc=NULL,impute=TRUE,select=TRUE,min.n=3,smooth=FALSE) { #function to score clusters according to the key matrix, correcting for item overlap
cl <- match.call() 
  
if(!inherits(stats, "statsBy")) stop("Please run statsBy first.  Make sure that you specified keys before stats.")
if(is.null(stats$r)) stop("I am sorry, but you need to run statsBy with the cors=TRUE option first")
MIN.n <- min.n
Item.smc <- item.smc  #we redefine later, so we need to keep this original value
bad.matrix <-0
na.matrix <- 0
r.list <- stats$r  #the call uses the output of statsBy 
Select <- select
keys.list<- keys #we are saving these for each pass	
ngroups <-length(r.list)	
 tol=sqrt(.Machine$double.eps)    #machine accuracy

var.names <- colnames(stats$r[[1]])

keys <- keys.list #these are the original keys   
select <- Select  #the original values
 bad <- FALSE
 if(is.list(keys) & (!is.data.frame(keys))) { if (select) {
    select <- selectFromKeyslist(var.names,keys)
 # select <- sub("-","",unlist(keys))   #added April 7, 2017
      select <- select[!duplicated(select)]
      }  else {select <- 1:length(var.names) }
      } 
 keys <- make.keys(stats$r[[1]][select,select],keys)  #added 9/9/16    (and then modified March 4, 2017  #move outside of loop
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix

#now, do repeated scoring, one pass for each group
result <- list(Call = cl)
cor.list <- list()
var.list <- list()
alpha.list <- list()
r.list <- list()
for(group in 1:ngroups) {

num.n <- stats$nWg[[group]] #this is the pairwise data 

#if(sum(!is.na(num.n)) <  10 || mean(num.n,na.rm=TRUE) < MIN.n){
if(mean(num.n,na.rm=TRUE) < MIN.n){
    cor.list[[group]] <- NA
    alpha.list[[group]] <- NA   #add an empty element
    r.list[[group]] <- NA
   next()
   }
r <- stats$r[[group]]
r[num.n < min.n] <- NA
  
 #if (!isCorrelation(r)) {r <- cor(r[,select],use="pairwise")} else 
r <- r[select,select]

  #only use correlations
 #if ((dim(r)[1] != dim(r)[2]) ) {r <- cor(r,use="pairwise")}
 #if(any(abs(r[!is.na(r)]) > 1)) warning("Something is seriously wrong with the correlation matrix, some correlations had absolute values > 1!  Please check your data.")
 if(any(is.na(r))) {
                SMC=FALSE
 #                warning("Missing values in the correlation matrix do not allow for SMC's to be found")
                 bad <- TRUE}

 if(SMC && is.null(Item.smc)) {item.smc <- smc(r)} else {
         diag(r) <- NA
          item.smc <- apply(r,1,function(x) max(abs(x),na.rm=TRUE))
          item.smc[is.infinite(item.smc) ] <- 1 
          diag(r) <- 1}
                                   
 if(all(item.smc ==1)) SMC <- FALSE

#item.smc <- rep(1,NCOL(r))  #why are we doing this?

 if(!bad) {covar <- t(keys) %*% r %*% keys} else  #matrix algebra is our friend 
     {
    covar<- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))*x,na.rm=TRUE))  #matrix multiplication without matrices!
   covar <- apply(keys,2,function(x)  colSums(apply(keys,2,function(x)  colMeans(r*x,na.rm=TRUE))*x,na.rm=TRUE)) *NROW(keys) 
 # we find the means of the non-zero elements by adjusting for the number of these elements
    covar <- score.na(keys,r,cor=FALSE)
  }
  
 
 
 n.keys <- ncol(keys)
 item.var <- item.smc
 raw.r  <- cov2cor(covar)
# key.var <- diag(t(keys) %*% keys)  #if correlations
 key.var <- diag(r) %*% abs(keys)   #if covariances or correlations
 key.n <-  diag(t(keys) %*% keys)
 var <- diag(covar)     #these are the scale variances unless we handled bad data
 key.smc <- as.vector( t(abs(keys)) %*% item.smc ) 
 #key.alpha <- ((var-key.var)/var)*(key.var/(key.var-1))
 key.alpha <- ((var-key.var)/var)*(key.n/(key.n-1))
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
        # item.cov <- apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))  #some correlations are NA which makes this value too low
          item.cov <- apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))    #we use the means to fix this problem 
         raw.cov <-  apply(keys,2,function(x) colSums(item.cov*x,na.rm=TRUE))
         #raw.cov <- covar
         item.cov <- t(item.cov)   #what is going on?
   }
 adj.cov <- raw.cov 


 #now adjust them
 
  med.r <- rep(NA, n.keys)
 for (i in 1:(n.keys)) {
    temp <- keys[,i][abs(keys[,i]) > 0]
   temp <- diag(temp,nrow=length(temp))
   small.r <- r[abs(keys[,i])>0,abs(keys[,i])>0]
   small.r <- temp %*% small.r %*% temp  #does not work with NA
    med.r[i]  <- median(small.r[lower.tri(small.r)],na.rm=TRUE)  
    for (j in 1:i) {
  
  #maybe don't adjust
   
if(av.r) { adj.cov[i,j] <- adj.cov[j,i]<- raw.cov[i,j] - sum(keys[,i] * keys[,j] ) + sum(keys[,i] * keys[,j] *  sqrt(key.av.r[i] * key.av.r[j]))
  } else {
     adj.cov[i,j] <- adj.cov[j,i] <- raw.cov[i,j] - sum(keys[,i] * keys[,j] )+ sum( keys[,i] * keys[,j] * sqrt(item.smc[i]* abs(keys[,i])*item.smc[j]*abs(keys[,j]) ))
 
 }
    } }

scale.var <- diag(raw.cov)

diag(adj.cov) <- diag(raw.cov)
 adj.r <- cov2cor(adj.cov)  #this is the overlap adjusted correlations

if(smooth) {
  if(any(is.na(adj.r))) {na.matrix <- na.matrix+1} else {
  ev <- eigen(adj.r)
  bad.matrix <- bad.matrix + any(ev$value<0)
  adj.r <- cor.smooth(adj.r)
  }
}

#adjust the item.cov for item overlap
#we do this by replacing the diagonal of the r matrix with the item.var (probably an smc, perhaps a maximum value)

diag(r) <- item.var
if(!bad) { item.cov <- t(keys) %*% r    #the normal case is to have all correlations
        } else {  
         #item.cov <- t(apply(keys,2,function(x) colSums(r*x,na.rm=TRUE)))  #some correlations are NA
              item.cov <- matrixMult.na(t(keys),r,scale=FALSE)
         }



 if(n.keys > 1) {
    item.cor <-   sqrt(diag(1/(key.lambda6*scale.var))) %*% (item.cov)  # %*% diag(1/sqrt(item.var))
    rownames(item.cor) <- colnames(keys)
    } else {
      item.cor <- r %*% keys /sqrt(key.lambda6*scale.var) }
   colnames(item.cor) <- colnames(r)
   item.cor <- t(item.cor)
   names(med.r) <- colnames(keys)


 

 adj.r.vect <- as.vector(adj.r[lower.tri(adj.r)])
 if (correct) {cluster.corrected <- correct.cor(adj.r,t(key.alpha))
 

r.list [[group]] <- adj.r.vect
 alpha.list[[group]] <-  list(alpha=key.alpha)
 cor.list[[group]] <- list(cor=adj.r)
 var.list [[group]] <- list(var=var/(key.n^2))
# names(result[group]) <- names(stats$r[group])
 }  #correct for attenuation
 else {

result[[group]] <- list(cor=adj.r,var=var, sd=sqrt(var),alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,item.cor=item.cor,med.r=med.r,Call=cl)}


}   #end of group loop



 names(alpha.list) <- names(cor.list) <- names(stats$r)
 if(!exists("adj.r.vect")) browser()                        #debugging option
# r.list <- r.list[!is.na(r.list)]  #this gets rid of the missing cases
 cor.mat <- matrix(unlist(r.list),ncol=length(adj.r.vect),byrow=TRUE)
 var.mat <- matrix(unlist(var.list),ncol=NCOL(keys),byrow=TRUE)
# rownames(cor.mat) <- names(cor.list[!is.na(cor.list)])

rownames(var.mat) <- names(cor.list[!is.na(cor.list)])
colnames(var.mat)<- colnames(keys)
 rownames(cor.mat) <- names(cor.list[!is.na(cor.list)])

rownames(cor.mat) <- names(cor.list)


 cnR <- abbreviate(colnames(keys),minlength=5) 
      k <- 1
      nvar <- NCOL(keys)
      temp.name <- rep(NA,nvar*(nvar-1)/2)  #strange, but seems necessary
     for(i in 1:(nvar-1)) {for (j in (i+1):nvar) {
     temp.name[k]  <- paste(cnR[i],cnR[j],sep="-") 
     # colnames(cor.mat)[k] <- paste(cnR[i],cnR[j],sep="-")
      k<- k +1 }}
      colnames(cor.mat)<- temp.name
 
 result <- list(alpha=alpha.list,cor = cor.list,cor.mat=cor.mat,bad.matrix=bad.matrix, na.matrix=na.matrix,var = var.mat)
 class(result) <- c ("psych", "scoreBy")
 return(result)}
 #modified 01/11/15 to find r if not a square matrix
#modifed 03/05/15 to do pseudo matrix multiplication in case of missing data 

select.from.list <- function(x){
ngroups <- length(x)
result <-list()
for (i in 1:ngroups) {
 result[[i]] <- x[[i]]$cor[lower.tri(x[[i]]$cor)]}
 return(result)
}