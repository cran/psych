"score.items"  <-
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE, impute="median",delete=TRUE,  min=NULL,max=NULL,digits=2) {
 message("score.items has been replaced by scoreItems, please change your call")
     scoreItems(keys=keys,items=items,totals=totals,ilabels=ilabels,missing=missing,impute=impute,delete=delete,min=min,max=max,digits=digits)
     }

"scoreItems"  <-
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE, impute="median",delete=TRUE,  min=NULL,max=NULL,digits=2) {
   cl <- match.call()
   raw.data <- TRUE
   keys <- as.matrix(keys)   #just in case they were not matrices to start with
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
     abskeys <- abs(keys)
     keynames <- colnames(keys)
   num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
  num.ob.item <- num.item   #will be adjusted in case of impute = FALSE
    if (!missing) items <-  na.omit(items) 
    n.subjects <- dim(items)[1]
     if ((dim(items)[1] == dim(items)[2])  &&  (((min(items,na.rm=TRUE) < -1) || (max(items,na.rm=TRUE) > 1)))) {warning("You have an equal number of rows and columns but do not seem to have  a correlation matrix.  I will treat this as a data matrix.")} # with the exception for the very unusual case of exactly as many items as cases reported by Jeromy Anglim 
    if ((dim(items)[1] == dim(items)[2])  &&  (!((min(items,na.rm=TRUE) < -1) || (max(items,na.rm=TRUE) > 1)))){ #this is the case of scoring correlation matrices instead of raw data  (checking for rare case as well)       
     raw.data <- FALSE
     n.subjects <- 0
     C <- as.matrix(items)
     cov.scales <- t(keys) %*% C %*% keys  #fast, but does not handle the problem of NA correlations
     cov.scales2 <- diag(t(abskeys) %*% C^2 %*% abskeys) # this is sum(C^2)  for finding ase
     response.freq <- NULL
           }  else {
   items <- as.matrix(items)
    
    response.freq <- response.frequencies(items)
    item.var <- apply(items,2,sd,na.rm=TRUE)
       bad <- which((item.var==0)|is.na(item.var))
       if((length(bad) > 0) && delete) {
       for (baddy in 1:length(bad)) {warning( "Item= ",colnames(items)[bad][baddy]  , " had no variance and was deleted from the data and the keys.")}
       items <- items[,-bad]
        keys <- as.matrix(keys[-bad,])
       
        n.items <- n.items - length(bad) 
        abskeys <- abs(keys)
        colnames(keys) <- keynames
      
        }
    item.means <- colMeans(items,na.rm=TRUE)
    if (is.null(min)) {min <- min(items,na.rm=TRUE)}
    if (is.null(max)) {max <- max(items,na.rm=TRUE)}
    # miss.rep <- rowSums(is.na(items))

     miss.rep <- (is.na(items) +0) %*% abs(keys)
    
   
    num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
    num.ob.item <- num.item   #will be adjusted in case of impute = FALSE
   if(impute !="none") {
        miss <- which(is.na(items),arr.ind=TRUE)
        if(impute=="mean") {
       		item.means <- colMeans(items,na.rm=TRUE)   #replace missing values with means
       		items[miss]<- item.means[miss[,2]]} else { 
       		item.med   <- apply(items,2,median,na.rm=TRUE) #replace missing with medians
        	items[miss]<- item.med[miss[,2]]}   #this only works if items is a matrix
         scores <- items %*%  keys  #this actually does all the work but doesn't handle missing values
          C <- cov(items,use="pairwise")
          cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
          cov.scales2 <- diag(t(abskeys) %*% C^2 %*% abskeys)   # sum(C^2)  for finding ase
        }  else { #handle the case of missing data without imputation
           scores <- matrix(NaN,ncol=n.keys,nrow=n.subjects)
           totals <- FALSE  #just in case it was not already false
           for (scale in 1:n.keys) {
           pos.item <- items[,which(keys[,scale] > 0)]
           neg.item <- items[,which(keys[,scale] < 0)]
           neg.item <- max + min - neg.item
           sub.item <- cbind(pos.item,neg.item)
           scores[,scale] <- rowMeans(sub.item,na.rm=TRUE)
           rs <- rowSums(!is.na(sub.item))
           num.ob.item[scale] <- mean(rs[rs>0])  #added Sept 15, 2011
          # num.ob.item[scale] <- mean(rowSums(!is.na(sub.item))) # dropped 
           		} # end of scale loop
       	
           # we now need to treat the data as if we had done correlations at input
            C <- cov(items,use="pairwise")
            cov.scales <- t(keys) %*% C %*% keys
            cov.scales2 <- diag(t(abskeys) %*% C^2 %*% abskeys)  # sum(C^2)  for finding ase
            raw.data <- FALSE
         }  #end of treating missing without imputation

       }
                   
    slabels <- colnames(keys)
    if (is.null(slabels)) {
    	if (totals) {slabels<- paste("S",1:n.keys,sep="")} else {
    	             slabels <- paste("A",1:n.keys,sep="")} }
   
    
   

    item.var <- diag(C)  #find the item variances
    
    var.scales <- diag(cov.scales)
    cor.scales <- cov2cor(cov.scales)    
    sum.item.var <- item.var %*% abskeys 
    sum.item.var2 <- item.var^2 %*% abskeys 
   
    
   #av.r <- (var.scales - sum.item.var)/(num.item*(num.item-1))  #actually, this the average covar
   alpha.scale <- (var.scales - sum.item.var)*num.item/((num.item-1)*var.scales)
   av.r <- alpha.scale/(num.item - alpha.scale*(num.item-1))  #alpha 1 = average r
   alpha.ob <- av.r * num.ob.item/(1+(num.ob.item-1)* av.r)
  
    colnames(alpha.scale) <- slabels
  alpha.scale[is.nan(alpha.scale)] <- 1
   
   #Find standard errors of alpha following Duhacheck and Iacobbci
   #Q = (2 * n^2/((n-1)^2*(sum(C)^3))) * (sum(C) * (tr(C^2) + (tr(C))^2) - 2*(tr(C) * sum(C^2)))
   #this works if we have the raw data
   Q = (2 * num.item^2/((num.item-1)^2*((var.scales)^3))) * (var.scales * (sum.item.var2 + sum.item.var^2) - 2* sum.item.var * cov.scales2)

   ase <- NULL  #to have something if we don't have raw data
   #now find the Guttman 6 * reliability estimate as well as the corrected item-whole correlations
   
   if(raw.data) { item.cor <- cor(items,scores)} else {if (n.keys >1) {
         item.cor <- C %*% keys %*% diag(1/sqrt(var.scales))/sqrt(item.var)} else {item.cor <- C %*% keys /sqrt(var.scales * item.var)}}
         colnames(item.cor) <- slabels
    c.smc <- smc(C,TRUE)
   
    diag(C) <- c.smc
    sum.smc <- c.smc %*% abskeys
    G6 <- (var.scales - sum.item.var + sum.smc)/var.scales
    corrected.var <- diag(t(keys) %*%  C %*% keys)
    if(n.keys>1) {
    item.rc <- (C %*% keys) %*% sqrt(diag(1/corrected.var))/sqrt(item.var)} else {
      item.rc <- C %*% keys /sqrt(corrected.var*item.var) }
    colnames(item.rc) <- slabels
   
  if(n.subjects >0) {ase <- sqrt(Q/ n.subjects )} else {ase=NULL}  #only meaningful if we have raw data
  if(is.null(ilabels)) {ilabels <- colnames(items) }
  if(is.null(ilabels)) {ilabels <-  paste("I",1:n.items,sep="")}
   rownames(item.rc) <- ilabels
  if(raw.data) {
     correction <- (colSums(abs(keys)-(keys))/2)*(max+min) #correct for flipping
     scores <- scores  + matrix(rep(correction,n.subjects),byrow=TRUE,nrow=n.subjects)
      
     if (!totals) {
        if(n.keys > 1) {scores <- scores %*% diag(1/num.item)   #find averages
                  }  else
                 {scores <- scores/num.item } }
            colnames(scores) <- slabels          
         } else {if (impute !="none") scores <- NULL}
    scale.cor <- correct.cor(cor.scales,t(alpha.scale))
  rownames(alpha.scale) <- "alpha"
  rownames(av.r) <- "average.r"
  rownames(G6) <- "Lambda.6"
  sn <-  av.r * num.item/(1-av.r)
  rownames(sn) <- "Signal/Noise"

   if (!raw.data) { 
     if(impute =="none") {
       rownames(alpha.ob) <- "alpha.observed"
       if(!is.null(scores)) colnames(scores) <- slabels #added Sept 23, 2013
       results <-list(scores=scores,missing = miss.rep,alpha=alpha.scale, av.r=av.r,sn=sn, n.items = num.item,  item.cor = item.cor,cor = cor.scales, corrected = scale.cor,G6=G6,item.corrected = item.rc,response.freq=response.freq,raw=FALSE,alpha.ob = alpha.ob,num.ob.item =num.ob.item,ase=ase,Call=cl)} else {
                            results <- list(alpha=alpha.scale, av.r=av.r,sn=sn, n.items = num.item,  item.cor = item.cor,cor = cor.scales ,corrected = scale.cor,G6=G6,item.corrected = item.rc ,response.freq =response.freq,raw=FALSE, ase=ase,Call=cl)}  } else {
   if(raw.data) {if (sum(miss.rep) > 0) {results <-list(scores=scores,missing = miss.rep,alpha=alpha.scale, av.r=av.r, sn=sn,n.items = num.item,  item.cor = item.cor,cor = cor.scales ,corrected = scale.cor,G6=G6,item.corrected = item.rc,response.freq=response.freq,raw=TRUE,ase=ase,Call=cl)} else{  
                                         results <- list(scores=scores,alpha=alpha.scale, av.r=av.r,sn=sn, n.items = num.item,  item.cor = item.cor, cor =cor.scales,corrected = scale.cor,G6=G6,item.corrected = item.rc ,response.freq=response.freq,raw=TRUE,ase=ase,Call=cl)} }
   }
   class(results) <- c("psych", "score.items")
    return(results)
 }
 #modified June 1 to add row names to items 
 #modified June 22 to add median imputation
 #modified August 8 to add colnames to scores
 #modified Sept 23, 2007 to allow for short output
 #modified December 10, 2007 to default to ilabels as colnames(items)
 #modified March, 2009 to better use the print.psych function
 #modified March 22, 2009 to add G6 and corrected for item overlap correlation
 #also allow for correlation/covariance matrix input
 #modified Sept 3, 2010 to include response frequencies
 #modified November 11, 2010 to allow for data with lots of missingness to be scored without imputing means or medians  
 #need to rethink the short option.  Why bother since summary and print don't show scores anyway
 #added missing score to count missing responses for each scale instead of just the overall.
 #Modified November 22, 2013 to add confidence intervals for alpha 