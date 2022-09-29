"score.items"  <-
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE, impute="median",delete=TRUE,  min=NULL,max=NULL,digits=2,select=TRUE) {
 message("score.items has been replaced by scoreItems, please change your call")
     scoreItems(keys=keys,items=items,totals=totals,ilabels=ilabels,missing=missing,impute=impute,delete=delete,min=min,max=max,digits=digits,select=select)
     }

"scoreItems"  <-
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE, impute="median",delete=TRUE,  min=NULL,max=NULL,digits=2,n.obs=NULL,select=TRUE) {
   cl <- match.call()
  
   ## First some housekeeping for various types of input
    raw.data <- TRUE
    
  # if(is.list(keys) & !is.data.frame(keys)) keys <- make.keys(items,keys)   #added 9/9/16  and then fixed March 4, following a suggestion by Jeromy Anglim
    if(is.null(colnames(items))  ) {select <- FALSE   #can not select items if they don't have colnames or if the keys don't have rownames
           colnames(items) <- paste0("V",1:NCOL(items)) }   #give the items some  names
    
     if(is.null(dim(keys)) &(is.null(names(keys)))) {keys <- as.matrix(keys)   #the case of unnamed key returned from alpha
      rownames(keys) <-colnames(items)
      colnames(keys) <- paste0("Scale1",1:NCOL(keys))} 
      
    if(is.matrix(keys) & is.null(colnames(keys))) { colnames(keys) <- paste0("Scale",1:NCOL(keys))
                        if(is.null(rownames(keys)))  rownames(keys) <- paste0("V",1:NCOL(items))}  #this handles input from GAabbreviate

    if (select) {if(is.list(keys) & (!is.data.frame(keys))) {
#  select <- sub("-","",unlist(keys))  #then, replaced with select option, Apri 8, 2017
  		select <- selectFromKeyslist(colnames(items),keys)
     	 select <- select[!duplicated(select)]
      }  else {keys <- keys2list(keys)
        select <- selectFromKeyslist(colnames(items),keys)
      	select <- select[!duplicated(select)]
         }   
# if (!isCorrelation(r)) {r <- cor(r[select],use="pairwise")} else {r <- r[select,select]}
#check for bad input   -- the Mollycoddle option 
if(any( !(select %in% colnames(items)) )) {
 cat("\nVariable names in keys are incorrectly specified. Offending items are ", select[which(!(select %in% colnames(items)))],"\n")
 stop("I am stopping because of improper input.   See above for a list of bad item(s). ")}
 keys <- make.keys(items[,select],keys)} else {select <- 1:ncol(items) }
   #modified once again April 6,2017 to allow for selecting items 
   if(is.list(keys)) keys <- make.keys(ncol(items),keys,item.labels=colnames(items))   #added 12/26/21 to handle select=FALSE 
   keys <- as.matrix(keys)   #just in case they were not matrices to start with
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
    abskeys <- abs(keys)
    keynames <- colnames(keys)
    num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
    num.ob.item <- num.item   #will be adjusted in case of impute = FALSE
    if (!missing) items <-  na.omit(items) 
    n.subjects <- dim(items)[1]
     if ((dim(items)[1] == dim(items)[2])  &&  !isCorrelation(items)) {warning("You have an equal number of rows and columns but do not seem to have  a correlation matrix.  I will treat this as a data matrix.")} # with the exception for the very unusual case of exactly as many items as cases reported by Jeromy Anglim 
    if ((dim(items)[1] == dim(items)[2])  &&  isCorrelation(items)){ #this is the case of scoring correlation matrices instead of raw data  (checking for rare case as well)       
     raw.data <- FALSE
     totals <- FALSE #because we don't have the raw data, totals would be meaningless
      items <- items[select,select]  #we have a correlation matrix, but we don't want all of it
     n.subjects <- 0
     C <- as.matrix(items)
     cov.scales <- t(keys) %*% C %*% keys  #fast, but does not handle the problem of NA correlations
     cov.scales2 <- diag(t(abskeys) %*% C^2 %*% abskeys) # this is sum(C^2)  for finding ase
     response.freq <- NULL
           }  else {
     items <- items[,select]  #select just the items that we want to score.  This is faster and robust to bad items
    #check to make sure all items are numeric  --  if not, convert them to numeric if possible, flagging the item that we have done so 
     if(!is.matrix(items)) {  #does not work for matrices
    for(i in 1:n.items) {   
        if(!is.numeric(items[[i]] ))  {
                               
                                  if(is.factor(unlist(items[[i]])) | is.character(unlist(items[[i]]))) {  items[[i]] <- as.numeric(items[[i]]) 
                                  
                                colnames(items)[i] <- paste0(colnames(items)[i],"*")
        
                          } else {items[[i]] <- NA} }
                         
               }
              } 
      

   items <- as.matrix(items)
    
    response.freq <- response.frequencies(items)
    item.var <- apply(items,2,sd,na.rm=TRUE)
       bad <- which((item.var==0)|is.na(item.var))   #is this a good idea? 
       #bad <- which(is.na(item.var)) 
 
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
          if(raw.data &&  totals == TRUE) warning("Specifying totals = TRUE without imputation can lead to serious problems.  Are you sure?")  #just in case it was not already false        #do we want to allow totals anyway?
           #we could try to parallelize this next loop
           for (scale in 1:n.keys) {
           	pos.item <- items[,which(keys[,scale] > 0)]
          	neg.item <- items[,which(keys[,scale] < 0)]
          	 neg.item <- max + min - neg.item
           	sub.item <- cbind(pos.item,neg.item)
           	if(!totals) {scores[,scale] <- rowMeans(sub.item,na.rm=TRUE)} else {scores[,scale] <- rowSums(sub.item,na.rm=TRUE)}
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

       }  #end of processing raw data 
                   
    slabels <- colnames(keys)
    if (is.null(slabels)) {
    	if (totals) {slabels<- paste("S",1:n.keys,sep="")} else {
    	             slabels <- paste("A",1:n.keys,sep="")} }
   
    
   

    item.var <- diag(C)  #find the item variances
    
    var.scales <- diag(cov.scales)
    cor.scales <- cov2cor(cov.scales)    
    sum.item.var <- item.var %*% abskeys 
    sum.item.var2 <- item.var^2 %*% abskeys 
    item.r <- cov2cor(C) #this does not handle pairwise complete correctly because cov and cor with pairwise work differently
  
   #find the median correlation within every scale 
   med.r <- rep(NA, n.keys)
   for(k in 1:n.keys) {
  
		temp <- keys[,k][abs(keys[,k]) > 0]
   		temp.keys <- temp
   		temp <- diag(temp,nrow=length(temp))
  # temp  <- diag(keys[,k ][abs(keys[,k])>0] )
	 small.r <- item.r[abs(keys[,k])>0,abs(keys[,k])>0,drop=FALSE]  #drop = FALSE for the case of 1 item scales
  # small.r <- temp %*% small.r %*% temp    #this flips, but does not work if any r is na
     small.r[temp.keys < 0,] <- -small.r[temp.keys< 0,,drop=FALSE]   #these two lines flip negatively keyed items even if some are NA (August 25, 2020)
      small.r[,temp.keys < 0] <- -small.r[,temp.keys< 0,drop=FALSE]
    	med.r[k]  <- median(small.r[lower.tri(small.r)],na.rm=TRUE)  
   } 
   
  names(med.r) <- slabels
    #but we want to do this for each scale
   
    
    
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
   
   if(raw.data) {
   if(length(bad) >0 ) {items <- items[,-bad]
    # item.var <-  item.var[-bad]
    # C <- C[-bad,-bad]
   } #this removes those items with no variance from the item statistics
    item.cor <- cor(items,scores)
    } else {if (n.keys >1) {
         item.cor <- C %*% keys %*% diag(1/sqrt(var.scales))/sqrt(item.var)} else {item.cor <- C %*% keys /sqrt(var.scales * item.var)}}
         colnames(item.cor) <- slabels
    c.smc <- smc(C,TRUE)
  
    diag(C) <- c.smc
    sum.smc <- c.smc %*% abskeys
    G6 <- (var.scales - sum.item.var + sum.smc)/var.scales
    corrected.cov <-  t(keys) %*%  C %*% keys
    corrected.var <- diag(corrected.cov)
    #corrected.var <- diag(t(keys) %*%  C %*% keys)
    if(n.keys>1) {
    item.rc <- (C %*% keys) %*% sqrt(diag(1/corrected.var))/sqrt(item.var)} else {
      item.rc <- C %*% keys /sqrt(corrected.var*item.var) }
    colnames(item.rc) <- slabels
    
    #put all of this much later 
    key.var <- diag(t(keys) %*% keys)
    # key.av.r <- key.alpha/(key.var - key.alpha*(key.var-1))  #alpha 1 = average r
     scale.size <- outer(key.var,key.var)
     MIMS <-  corrected.cov/scale.size
     diag(MIMS)<- av.r   
    
  if(n.subjects > 0) {ase <- sqrt(Q/ n.subjects )} else {if(!is.null(n.obs)) {ase <- sqrt(Q/ n.obs )} else {ase=NULL}}  #only meaningful if we have raw data
  if(is.null(ilabels)) {ilabels <- colnames(items) }
  if(is.null(ilabels)) {ilabels <-  paste("I",1:n.items,sep="")}
   #if(length(bad) >0 ) {rownames(item.rc) <- ilabels[-bad]} else {rownames(item.rc) <- ilabels}
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
 # rownames(med.r) <- "median.r"
  rownames(G6) <- "Lambda.6"
  sn <-  av.r * num.item/(1-av.r)
  rownames(sn) <- "Signal/Noise"
  
  
  #find the Multi-Item Multi Trait item x scale correlations   # added 9/5/22
 MIMT <- matrix(NA,n.keys,n.keys)

 for (i in 1:(n.keys)) {
    temp <- keys[,i][abs(keys[,i]) > 0]
    flip.item <- temp * item.cor[names(temp),,drop=FALSE]
   if(length(names(temp)) > 1) { MIMT[i,] <- colMeans(flip.item[names(temp),,drop=FALSE])} else {MIMT[i,] <- flip.item}
    }
  colnames(MIMT) <- rownames(MIMT) <- colnames(keys)
 keys.list <- keys2list(keys)  #put them into a list to save them
   if (!raw.data) { 
     if(impute =="none") {
       #rownames(alpha.ob) <- "alpha.observed"
       if(!is.null(scores)) colnames(scores) <- slabels #added Sept 23, 2013
       results <- list(scores=scores,missing = miss.rep,alpha=alpha.scale, av.r=av.r,sn=sn, n.items = num.item,  item.cor = item.cor,cor = cor.scales, corrected = scale.cor,G6=G6,item.corrected = item.rc,response.freq=response.freq,raw=FALSE,alpha.ob = alpha.ob,num.ob.item =num.ob.item,ase=ase,med.r=med.r,keys=keys.list,MIMS=MIMS,MIMT=MIMT,Call=cl)} else {
         results <- list(alpha=alpha.scale, av.r=av.r,sn=sn, n.items = num.item,  item.cor = item.cor,cor = cor.scales ,corrected = scale.cor,G6=G6,item.corrected = item.rc ,response.freq =response.freq,raw=FALSE, ase=ase,med.r=med.r,keys=keys.list,MIMS=MIMS,MIMT=MIMT, all=cl)}  } else {
   if(raw.data) {if (sum(miss.rep) > 0) {results <-list(scores=scores,missing = miss.rep,alpha=alpha.scale, av.r=av.r, sn=sn,n.items = num.item,  item.cor = item.cor,cor = cor.scales ,corrected = scale.cor,G6=G6,item.corrected = item.rc,response.freq=response.freq,raw=TRUE,ase=ase,med.r=med.r,keys=keys.list,MIMS=MIMS,MIMT=MIMT,Call=cl)} else{  
                                         results <- list(scores=scores,alpha=alpha.scale, av.r=av.r,sn=sn, n.items = num.item,  item.cor = item.cor, cor =cor.scales,corrected = scale.cor,G6=G6,item.corrected = item.rc ,response.freq=response.freq,raw=TRUE,ase=ase,med.r=med.r,keys=keys.list,MIMS=MIMS,MIMT=MIMT,Call=cl)} }
   }
   class(results) <- c("psych", "score.items")
    return(results)
 }
 #modified June 1 to add row names to items 
 #modified June 22 to add median imputation
 #modified August 8 to add colnames to scores
 #modified Sept 23, 2007 to allow for short output
 #modified December 10, 2007 to default to ilabels as colnames(items)
 #modified March, 2009 to better#Find scores from correlations for each subject
#impute missing values?
#The correlations for each subject are found by statsBy
#and are stored as a list of matrices in the r element

"scoreBy" <-
function(keys,stats,correct=TRUE,SMC=TRUE,av.r=TRUE,item.smc=NULL,impute=TRUE,select=TRUE) { #function to score clusters according to the key matrix, correcting for item overlap
cl <- match.call() 
if(!inherits(stats, "statsBy")) stop("Please run statsBy first")
if(is.null(stats$r)) stop("I am sorry, but you need to run statsBy with the cors=TRUE option first")
MIN.n <- 3
r.list <- stats$r  #the call uses the output of statsBy 
Select <- select
keys.list<- keys #we are saving these for each pass	
ngroups <-length(r.list)	
 tol=sqrt(.Machine$double.eps)    #machine accuracy

#now, do repeated scoring, one pass for each group
result <- list(Call = cl)
cor.list <- list()
alpha.list <- list()
r.list <- list()
for(group in 1:ngroups) {

num.n <- stats$n[group,]

#if(sum(!is.na(num.n)) <  10 || mean(num.n,na.rm=TRUE) < MIN.n){
if(mean(num.n,na.rm=TRUE) < MIN.n){
    cor.list[[group]] <- NA
    alpha.list[[group]] <- NA   #add an empty element
    r.list[[group]] <- NA
   next()
   }
r <- stats$r[[group]]


keys <- keys.list #these are the original keys
select <- Select  #the original values
 bad <- FALSE
 if(is.list(keys) & (!is.data.frame(keys))) { if (select) {
    select <- selectFromKeyslist(colnames(r),keys)
 # select <- sub("-","",unlist(keys))   #added April 7, 2017
      select <- select[!duplicated(select)]
      }  else {select <- 1:ncol(r) }   
 if (!isCorrelation(r)) {r <- cor(r[,select],use="pairwise")} else {r <- r[select,select]}
 keys <- make.keys(r,keys)}  #added 9/9/16    (and then modified March 4, 2017
 if(!is.matrix(keys)) keys <- as.matrix(keys)  #keys are sometimes a data frame - must be a matrix
  #only use correlations
 #if ((dim(r)[1] != dim(r)[2]) ) {r <- cor(r,use="pairwise")}
 #if(any(abs(r[!is.na(r)]) > 1)) warning("Something is seriously wrong with the correlation matrix, some correlations had absolute values > 1!  Please check your data.")
 if(any(is.na(r))) {
                SMC=FALSE
 #                warning("Missing values in the correlation matrix do not allow for SMC's to be found")
                 bad <- TRUE}

 if(SMC && is.null(item.smc)) {item.smc <- smc(r)} else {
         diag(r) <- NA
         item.smc <- apply(r,1,function(x) max(abs(x),na.rm=TRUE))
         item.smc[is.infinite(item.smc) ] <- 1 
         diag(r) <- 1}
                                   
 if(all(item.smc ==1)) SMC <- FALSE
 if(!bad) {covar <- t(keys) %*% r %*% keys} else  #matrix algebra is our friend 
     {covar<- apply(keys,2,function(x) colSums(apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))*x,na.rm=TRUE))  #matrix multiplication without matrices!
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
         item.cov <- apply(keys,2,function(x) colSums(r*x,na.rm=TRUE))  #some correlations are NA
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
   small.r <- temp %*% small.r %*% temp
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

#adjust the item.cov for item overlap
#we do this by replacing the diagonal of the r matrix with the item.var (probably an smc, perhaps a maximum value)

diag(r) <- item.var
if(!bad) { item.cov <- t(keys) %*% r    #the normal case is to have all correlations
        } else {  
         item.cov <- t(apply(keys,2,function(x) colSums(r*x,na.rm=TRUE)))  #some correlations are NA
         }



 if(n.keys > 1) {
    item.cor <-   sqrt(diag(1/(key.lambda6*scale.var))) %*% (item.cov)  # %*% diag(1/sqrt(item.var))
    rownames(item.cor) <- colnames(keys)
    } else {
      item.cor <- r %*% keys /sqrt(key.lambda6*scale.var) }
   colnames(item.cor) <- colnames(r)
   item.cor <- t(item.cor)
   names(med.r) <- colnames(keys)


 
 
 
 if (correct) {cluster.corrected <- correct.cor(adj.r,t(key.alpha))
 
adj.r.vect <- as.vector(adj.r[lower.tri(adj.r)])
r.list [[group]] <- adj.r.vect
 alpha.list[[group]] <-  list(alpha=key.alpha)
 cor.list[[group]] <- list(cor=adj.r)
# names(result[group]) <- names(stats$r[group])
 }  #correct for attenuation
 else {
result[[group]] <- list(cor=adj.r,sd=sqrt(var),alpha=key.alpha,av.r = key.av.r,size=key.var,sn=sn,G6 =key.lambda6,item.cor=item.cor,med.r=med.r,Call=cl)}


}   #end of group loop
 names(alpha.list) <- names(cor.list) <- names(stats$r)

 r.list <- r.list[!is.na(r.list)]  #this gets rid of the missing cases
 cor.mat <- matrix(unlist(r.list),ncol=length(adj.r.vect),byrow=TRUE)

# rownames(cor.mat) <- names(cor.list[!is.na(cor.list)])
 
 cnR <- abbreviate(colnames(keys),minlength=5) 
      k <- 1
      nvar <- NCOL(keys)
      temp.name <- rep(NA,nvar*(nvar-1)/2)  #strange, but seems necessary
     for(i in 1:(nvar-1)) {for (j in (i+1):nvar) {
     temp.name[k]  <- paste(cnR[i],cnR[j],sep="-") 
     # colnames(cor.mat)[k] <- paste(cnR[i],cnR[j],sep="-")
      k<- k +1 }}
      colnames(cor.mat)<- temp.name
 
 result <- list(alpha=alpha.list,cor = cor.list,cor.mat=cor.mat)
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
 #use the print.psych function
 #modified March 22, 2009 to add G6 and corrected for item overlap correlation
 #also allow for correlation/covariance matrix input
 #modified Sept 3, 2010 to include response frequencies
 #modified November 11, 2010 to allow for data with lots of missingness to be scored without imputing means or medians  
 #need to rethink the short option.  Why bother since summary and print don't show scores anyway
 #added missing score to count missing responses for each scale instead of just the overall.
 #Modified November 22, 2013 to add confidence intervals for alpha 
 #modified Sept 9, 2016 to add the keys.list option for the scoring keys
 #modified April 22, 2018 to include median within scale correlation
 #Modified August 25, 2020 to calculate median correlation in the case of variables with no variance
 
 "selectFromKeyslist" <- function(itemname,keys) {nkey <- length(keys)
 select <- NULL
   for (key in 1:nkey) {
   if(is.null(keys[[key]])) {select <- NULL} else {
      if(is.numeric(keys[[key]])) {select <- c(select,itemname[abs(unlist(keys[[key]]))]) } else {select <- c(select,sub("-", "", unlist(keys[[key]]))) }    
    }}
    return(select)}
      
  