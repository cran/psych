 #Pearson or polychoric correlations with confidence intervals

 "cor.ci" <- 
function(x, keys = NULL, n.iter = 100, p = 0.05, overlap=FALSE, poly = FALSE, method = "pearson",plot=TRUE,minlength=5,n=NULL,...) {
corCi(x=x, keys = keys, n.iter = n.iter, p = p, overlap=overlap, poly = poly, method = method,plot=plot,minlength=minlength,n=n,...) }

 "corCi" <- 
function(x,keys = NULL, n.iter = 100, p = 0.05, overlap=FALSE, poly = FALSE, method = "pearson",plot=TRUE,minlength=5,n = NULL,...) {

 cl <- match.call()
 n.obs <- dim(x)[1]
 
 if(!isCorrelation(x)) {#the normal case is to have data and find the correlations and then bootstrap them
 #added the direct from correlation matrix option, August 17, 2019 since I was finding them for statsBy
 
 if(is.null(keys) && overlap) overlap <- FALSE  #can not correct for overlap with just items
if(poly) {  #find polychoric or tetrachoric correlations if desired
  ncat <- 8
  nvar <- dim(x)[2]
  tx <- table(as.matrix(x))
	
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x)
	    typ = "tet"} else {  #should we do mixed correlations?
	    tab <- apply(x,2,function(x) table(x))
       if(is.list(tab)) {len <- lapply(tab,function(x) length(x))} else {len <- dim(tab)[1] }

dvars <- subset(1:nvar,len==2)   #find the dichotomous variables
pvars <- subset(1:nvar,((len > 2) & (len <= ncat)))  #find the polytomous variables
cvars <- subset(1:nvar,(len > ncat))  #find the continuous variables (more than ncat levels)
if(length(pvars)==ncol(x)) {tet <- polychoric(x)
	    typ = "poly"} else {tet <- mixedCor(x)
	    typ="mixed" }	    
	  } 
   
    Rho <- tet$rho    #Rho is the big correlation of all of items 
 } else { Rho <- cor(x,use="pairwise",method=method)    #the normal Pearson correlations
         } 
  #now, if there are keys, find the correlations of the scales       
   if(!is.null(keys)) {bad <- FALSE
   if(!is.matrix(keys)) keys <- make.keys(x,keys)  #handles the new normal way of just passing a keys list
         if(any(is.na(Rho))) {warning(sum(is.na(Rho)), " of the item correlations are NA and thus finding scales that include those items will not work.\n We will try to do it for those  scales which do not include those items.
         \n Proceed with caution. ") 
         bad <- TRUE
         rho <- apply(keys,2,function(x) colMeans(apply(keys,2,function(x) colMeans(Rho*x,na.rm=TRUE))*x,na.rm=TRUE))  #matrix multiplication without matrices!  
         #switched to using colMeans instead of colSums, recognizing the problem of different number of items being dropped.
         } else {
          rho <- t(keys) %*% Rho %*% keys} }  else {rho <- Rho} #find the covariances between the scales

 #

 ##correct for overlap if necessary on the original data

 if(overlap) { key.var <- diag(t(keys) %*% keys)
  var <- diag(rho)    #these are the scale variances
   n.keys <- ncol(keys)
   key.av.r <- (var - key.var)/(key.var * (key.var-1))
   item.cov <- t(keys) %*% Rho #this will blow up if there are bad data
   raw.cov <- item.cov %*% keys
   adj.cov <- raw.cov
    for (i in 1:(n.keys)) {
    for (j in 1:i) {
         adj.cov[i,j] <- adj.cov[j,i]<- raw.cov[i,j] - sum(keys[,i] * keys[,j] ) + sum(keys[,i] * keys[,j] *  sqrt(key.av.r[i] * key.av.r[j]))
                    }
                           }
       diag(adj.cov) <- diag(raw.cov)
       rho <- cov2cor(adj.cov)
     }
  rho <- cov2cor(rho)     #scale covariances to correlations  
  nvar <- dim(rho)[2]
 
  
 if(n.iter > 1) {
 	replicates <- list()
 	rep.rots <- list()
 	##now replicate it to get confidence intervals
   #replicates <- lapply(1:n.iter,function(XX) {
   
 		replicates <- mclapply(1:n.iter,function(XX) {
		 xs <- x[sample(n.obs,n.obs,replace=TRUE),]
 		{if(poly) {
 		 if(typ!= "tet") {tets <- mixedCor(xs)} else {tets <- tetrachoric(xs)}

  	R <- tets$rho} else {R <- cor(xs,use="pairwise",method=method)}  #R is the big correlation matrix
 
 if(!is.null(keys)) { if (bad) {covariances <- apply(keys,2,function(x) colMeans(apply(keys,2,function(x) colMeans(R*x,na.rm=TRUE))*x,na.rm=TRUE))  #matrix multiplication without matrices!
  } else {
 covariances <- t(keys) %*% R %*% keys}
               r <- cov2cor(covariances) } else {r <- R}					
  #correct for overlap if this is requested
  if(overlap) { 
  var <- diag(covariances) 
   item.cov <- t(keys) %*% R 
   raw.cov <- item.cov %*% keys
   adj.cov <- raw.cov
   key.av.r <- (var - key.var)/(key.var * (key.var-1))
    for (i in 1:(n.keys)) {
    for (j in 1:i) {
         adj.cov[i,j] <- adj.cov[j,i]<- raw.cov[i,j] - sum(keys[,i] * keys[,j] ) + sum(keys[,i] * keys[,j] *  sqrt(key.av.r[i] * key.av.r[j]))
                    }
                           }
       diag(adj.cov) <- diag(raw.cov)
       r <- cov2cor(adj.cov) #fixed 03/12/14
  }    
 rep.rots <- r[lower.tri(r)]
 }
 } 
 )
 
 }

   replicates <- matrix(unlist(replicates),ncol=nvar*(nvar-1)/2,byrow=TRUE) 
       z.rot <- fisherz(replicates)
     means.rot <- colMeans(z.rot,na.rm=TRUE)
      sds.rot <- apply(z.rot,2,sd, na.rm=TRUE)
      sds.rot <- fisherz2r(sds.rot)
      ci.rot.lower <- means.rot + qnorm(p/2) * sds.rot  #This is the normal value of the observed distribution
      ci.rot.upper <- means.rot + qnorm(1-p/2) * sds.rot
      means.rot <- fisherz2r(means.rot)
      ci.rot.lower <- fisherz2r(ci.rot.lower)
      ci.rot.upper <- fisherz2r(ci.rot.upper)
      low.e <- apply(replicates,2,quantile, p/2,na.rm=TRUE)
      up.e  <- apply(replicates, 2, quantile, 1-p/2,na.rm=TRUE)
      tci <- abs(means.rot)/sds.rot
      ptci <- 1- pnorm(tci)    #subtract from 1 (added 11/14/20)
      ci.rot <- data.frame(lower=ci.rot.lower,low.e=low.e,upper=ci.rot.upper,up.e=up.e,p =2*(ptci))   #dropped the 1-ptci 11/14/20
      cnR <- abbreviate(colnames(rho),minlength=minlength) 
      k <- 1
     for(i in 1:(nvar-1)) {for (j in (i+1):nvar) {
      rownames(ci.rot)[k] <- paste(cnR[i],cnR[j],sep="-")
      k<- k +1 }}
      
     


results <- list(rho=rho, means=means.rot,sds=sds.rot,tci=tci,ptci=ptci,ci=ci.rot,Call= cl,replicates=replicates)
#if(plot) {cor.plot.upperLowerCi(results,numbers=TRUE,cuts=c(.001,.01,.05),...) }  #automatically plot the results
if(plot) {cor.plot(rho,numbers=TRUE,cuts=c(.001,.01,.05),pval =  2*(1-ptci),...) }
class(results) <- c("psych","cor.ci")
    
return(results)  } else {#we have been given correlations, just find the cis.
if(is.null(n)) {warning("\nFinding confidence intervals from a correlation matrix, but  n is not specified, arbitrarily set to 1000")
   n <- 1000}
results <- cor.Ci(x,n=n, alpha=p, minlength=minlength)
results$ci <- results$r.ci
results$r <- x
class(results) <- cs(psych, corCi)
return(results)
}

 }
 #written Sept 20, 2013 
 #adapted from fa.poly
 #modified May 1, 2014 to scale by pvals
 #modified August 24, 2017 to include Bonferoni corrections from cor.test
 #modified January 22, 2023 to find correlations first if given raw data
 "corPlotUpperLowerCi" <-  "cor.plot.upperLowerCi" <- 
function(R,numbers=TRUE,cuts=c(.001,.01,.05),select=NULL,main="Upper and lower confidence intervals of correlations",adjust=FALSE,...) {
cor.cip <- NULL
names <- cs(corCi, cor.cip, corr.test ,cor.ci)
    value <- inherits(R,names,which=TRUE)  # value <- class(x)[2]
if(max(value)==0) {#find the correlations first if it is not a correlation matrix
                   if(isCorrelation(R)) {R <- corCi(R)} else {if(NROW(R)!= NCOL(R)) {R <-corr.test(R)}
                   # }  else { browser()
  # stop("I am stopping because the input is incorrect \ndata must either be a raw data matrix or the output of corr.test.")} 
               }}
if(adjust) {lower <- R$ci.adj$lower.adj
   upper <- R$ci.adj$upper.adj} else {
   lower <- R$ci$lower
upper <- R$ci$upper}
temp <- lower 
if(is.null(R$r)) {cn = colnames(R$rho)
         rl <- R$rho[lower.tri(R$rho)]} else {
         cn = colnames(R$r)
         rl <-  R$r[lower.tri(R$r)]} #is input from cor.ci or corr.test
lower[rl < 0 ] <- upper[rl < 0]
upper[rl < 0] <- temp[rl < 0]
m <- length(lower)
n <- floor((sqrt(1 + 8 * m) +1)/2) 
    X <- diag(n)
    X[lower.tri(X)] <- upper
    X <- t(X) 
    X[lower.tri(X)] <- lower
    diag(X) <- 1 
colnames(X) <- rownames(X) <- cn
if(is.null(R$ptci))  {pval <- R$p} else {pval = 2*(1-R$ptci)}

corPlot(X,numbers=numbers,pval=pval,cuts=cuts,select=select,main=main,...)  
class(X) <-  c("psych","cor.cip")
colnames(X) <- abbreviate(rownames(X,4))
invisible(X) }
 

 