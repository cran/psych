#added  just correlate with criteria to speed it up (June 23, 2017)
"bestItems" <- 
function(x,criteria=1,cut=.3, abs=TRUE, dictionary=NULL,check=FALSE,digits=2) {

if(check) {item.var <- apply(x,2,sd,na.rm=TRUE)  #added check 10/14/17
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0) ) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
             }
             }
if(NROW(criteria)> 1)  {x <- cbind(x,criteria)   #to allow for a separate object
    criteria <- "criteria" }
    
if(!isCorrelation(x))  { x <- cor(x,x[criteria],use="pairwise")} #the normal case --convert to correlation if necessary
  
if(abs) {ord <- order(abs(x[,criteria]),decreasing=TRUE)
  value <- x[ord,criteria,drop=FALSE]
  count <- sum(abs(value) > cut,na.rm=TRUE)
  value <- value[1:count,,drop=FALSE]
  } else {ord <- order(x[,criteria],decreasing=TRUE)
  value <- x[ord,criteria]
  value <- value[value,criteria > cut] }
value <- round(data.frame(value),digits)
if((!is.null(dictionary)) && !is.factor(dictionary)) {temp <- lookup(rownames(value),dictionary)
   value <- merge(value,temp,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(value) <- value[,"Row.names"]
   value <- value[-1]
  if(abs) {ord <- order(abs(value[,criteria]),decreasing=TRUE) } else {ord <- order(value[,criteria],decreasing=TRUE)}
   value <- value[ord,] 
   }
return(value)
}
  
  
  #lookup which x's are found in y[c1],return matches for y[]
"lookup" <- 
function(x,y,criteria=NULL) {
if (is.null(criteria)) {temp <- match(x,rownames(y))} else {
     temp <- match(x,y[,criteria])}
 y <- (y[temp[!is.na(temp)],,drop=FALSE])
return(y)}
 
  

"bestScales" <- 
 function(x,criteria,cut=.1,n.item =10,overlap=FALSE,dictionary=NULL,check=FALSE,impute="none", n.iter =1,frac=.9,digits=2) {
 cl <- match.call()
  first <- TRUE
if(check) {item.var <- apply(x,2,sd,na.rm=TRUE)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0) ) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
             }
             }
  #first, define function to be parallelized
   #mcmapply for parallel, mapply for debugging 
   ##
short <- function(i,x,n.obs,criteria,cut,n.item,impute,digits,dictionary,frac) {
   
	if(n.iter > 1) {
	   if(!isCorrelation(x)) { ss <- (1:n.obs) 
	    ss <- sample(ss,n.obs,FALSE) 
	    ss <- ss[1:(n.obs*frac)] 
	   scores   <-  bScales(x[ss,],criteria=criteria,cut=cut,
	             n.item =n.item,overlap=overlap,dictionary=dictionary,impute=impute,digits=digits)  
	             } else {message("iterative solutions not possible for correlation matrices")
	            n.iter <- 1 
	            }} else { # a correlation matrix or n.iter = 1
	scores   <-  bScales(x,criteria=criteria,cut=cut,
	             n.item =n.item,overlap=overlap,dictionary=dictionary,impute=impute,digits=digits)
	             }
	
	 key.list <- keys2list(scores$key)
	 if((n.iter > 1) & !first) {
  		cross <- scoreFast(key.list,x[-ss,],impute=impute,min=1,max=6) 
 	 validity <- diag(cor(cross,x[-ss,criteria],use="pairwise"))
  	result <- list(validity=c(scores$r,validity),keys=key.list,R = scores$R)
  } else {result <- scores
          result$key.list <- key.list}
          class(result) <- cbind("psych","bestScales")
 return(result) }
 ###
 
#begin the main function 
#if criteria is a separate data frame, combine x and criteria
#there are several alternative forms for criteria
#it is either a column name of x, or it is a separate data.frame/matrix
if(!is.null(dim(criteria))| (length(criteria) == NROW(x)))  { x <- cbind(x,criteria)
    if(NCOL(criteria) > 1 ){criteria <- colnames(criteria) } else {criteria <- "criteria"}
    }
 
 n.obs <- nrow(x)
 if((n.iter ==1)| first ) {   #don't bother to parallelize, just do one trial
   first.result <- short(1,x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=1)
   first <- FALSE
   result <- first.result
   }
     #the case for n.iter > 1.  We want to parallelize this because we are working pretty hard
 if(n.iter > 1) { 
result <- list()

#make mcmapply when not debgging
result <- mcmapply(short,c(1:n.iter),MoreArgs=list(x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=frac))
  
 #save the keys and the summary 
  validity <- list()
  #validity is a list of 3 elements repeated n.iter times
  #first are the validities
  #then are the keys
  #then are the item by criteria correlations

  keys <- list()
  R.list <- list()
  
  for(i in (1:n.iter)) { validity[[i]] <- result[[i*3 -2]]
    keys[[i]] <- result[[i*3-1]]
    R.list[[i]] <- result[[i*3]]
    }
     

   replicated.items <- bestReplicatedItems(keys)

   items <- list()
    item.mean <- list()
   for(j in 1:length(criteria)) {
     #first, find  the means and standard deviations for each selected item
      if(length(criteria) > 1 ) {for (i in 1:n.iter) { item.mean[[i]] <-  R.list[[i]][names(replicated.items[[j]]),criteria[j]] }
      } else { for (i in 1:n.iter) {item.mean[[i]] <- R.list[[i]][names(replicated.items[[j]])] } }
     item.m <- matrix(unlist(item.mean),nrow=n.iter,byrow=TRUE)
   
   
     colnames(item.m) <- names(replicated.items[[j]])
      means = colMeans(item.m,na.rm=TRUE)
     sds <- apply(item.m,2,sd,na.rm=TRUE)  
	items [[criteria[j] ]] <-  cbind(replicated.items[[j]],mean.r=means,sd.r = sds,dictionary[names(replicated.items[[j]]),])
	items[[criteria[j]]] <- dfOrder(items [[criteria[j] ]],"-mean.r",absolute=TRUE) #sort on the mean.r column

    }
    				

   result.df <- data.frame(matrix(unlist(validity),ncol=2*length(criteria),byrow=TRUE))
  colnames(result.df) <-c(paste("derivation",criteria),paste("validation",criteria)) 
  
  result <- list(validity = result.df,items=items,replicated.items =replicated.items,keys = keys,Call=cl)
  result$first.result  <- first.result
  result$means = colMeans(result.df)
  result$sd <- apply(result.df,2,sd)
  ncriteria <- length(criteria)
  result$summary <- data.frame(derivation.mean = result$means[1:ncriteria],derivation.sd = result$sd[1:ncriteria],validation.m=result$mean[(ncriteria+1):(2*ncriteria)],
            validation.sd =result$sd[(ncriteria+1):(2*ncriteria)] )
    rownames(result$summary) <- criteria
  
   }
  
  class(result) <- c("psych","bestScales")
  return(result)
}

"bScales" <- 
function(x,criteria,cut=.1,n.item =10, overlap=FALSE,dictionary=NULL,impute="median",digits=2) {

 #created 20/2/14
 #find the scales based upon the items that most correlate with a criteria
 #pure dust bowl empiricism
 #modified 13/3/15 to handle the problem of missing item labels
 #Completely redone June 2017 to allow for raw data and bootstrapping
 ##


#first, declare a function to identify the bad items and drop them from the keys
 findBad <- function(key,r) { 
	ss <- abs(key) > 0 
	rss <- r[ss,ss] 
	if(any(is.na(rss))){ #some of these are bad
	n.bad <-  apply(rss,1,function(x) sum(is.na(x)))
	key[names(which.max(n.bad))] <- 0
	findBad(key,r)}
	return(key)
	}

short <- function(key,r) { 
	 kn <- names(key[abs(key[,1]) > 0,1])
	if(is.null(kn)) kn <- names(which(abs(key[,1]) > 0))
	 cn <- colnames(key)
 	ord <- order(abs(r[kn,cn]),decreasing=TRUE)
 	kn <- kn[ord]
	 result <- r[kn,cn,drop=FALSE]
	 return(result)
	}

#begin the main function
 ##Basically two cases:
 #a correlation matrix is provided and we do basic matrix algebra
 #raw data is provided (getting ready to do bootstrapping) and we find just the necessary correlations
 
 nvar <- ncol(x)
 if(isCorrelation(x)) {r <- x      #  case 1
    raw <- FALSE} else {  #case 2
    y <- x[,criteria]
    r <- cor(x,y,use="pairwise")
    colnames(r) <- criteria
     x <- as.matrix(x)
     raw <- TRUE
     n.obs <- NROW(x)}
    #don't actually need to have  a square matrix
 ny <- length(criteria)
 nc <- length(cut)
 ni <- length(n.item)   #number of items per scale to find
 ord.name <- NULL
if(length(cut) == 1)  cut <- rep(cut,ny)
if(length(n.item) == 1) n.item <- rep(n.item,ny)

#this next part just finds the cut values to use
 if(!overlap)  {r[criteria,criteria] <- 0} else {for(i in 1:ny) r[criteria[i],criteria[i]] <- 0}
 if(ny > 1 ) {ord <- apply(abs(r[,criteria]),2,order,decreasing=TRUE) 
     for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[ord[n.item[i],i],criteria[i]])) 
    # ord.name <- c(ord.name, rownames(r)[ord[1:n.item[i],i]] )     #this is wrong and not needed
    }
     } else {
         ord <- order(abs(r[,criteria]),decreasing=TRUE)
         for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[ord[n.item[i]+1],criteria])) }
        }
#    cut has been adjusted

 key <- matrix(0,ncol=ny,nrow=nvar)
 key[t(t(r[,criteria]) >= cut)] <- 1
 key[t(t(r[,criteria]) <= -cut)]<- -1
 rownames(key)  <- rownames(r)
 colnames(key)  <- criteria

 k <- key  #this just gets it to be a matrix of the right size and names

 #colnames(key) <- paste(criteria,"S",sep=".")
 colnames(key) <- criteria
 #now, drop those items from the keys that are not used
 used <- rowSums(abs(key))
 key <- key[used >0,,drop=FALSE]  
 x <- x[,used >0,drop=FALSE]
 
#now, if we have raw data, find the correlation of the composite scale with the criteria
#if we have raw data, then we find the scales from the data 
if(raw)  { #case 2
#score <- matrix(NA,ncol=ny,nrow=nrow(x))
#for (i in (1:ny)) {
#   score[,i] <- rowSums(t((key[,i])* t(x)),na.rm=TRUE)
 #     }
 score <- scoreFast(key,x,impute=impute,min=1,max=6)
   R <- diag(cor(score,y,use="pairwise"))  #the validities
re <- r[,criteria]
ni <- colSums(abs(key))




} else {  #case 1 (from a correlation matrix)
  score <-NULL
  r <- r[,used > 0,drop=FALSE]
	if(any(is.na(r))) {#Are there any bad values
 		 for(i in 1:ny) {#key[,i] <- findBad(key[,i],r)  #Drop the bad items from any scoring key
 		 k[,i] <- colSums(t((key[,i]) * t(r)),na.rm=TRUE)}    #replace matrix addition with a colSums
 		 k <- t(k)
	} else {#otherwise, don't bother
    C <-  t(t(key) %*% t(r[criteria,,drop=FALSE]))  #criterion covariance
    V <-  t(key) %*%  r[ used > 0,] %*% key   #predictor variance
    
#	k  <- t(t(key) %*% t(r[criteria,,drop=FALSE]))  #we can do the matrix multiply because there are no bad data         
 	}
#	V <- t(k) %*% key   #this is the covariance of the criteria with criteria
#	C <- k[criteria,] 
 if(ny < 2) {re <- r[criteria,] 
           R <- C/sqrt(V)} else {
   R <- diag(C/sqrt(V))
	#re <- diag(k[criteria,])/sqrt(diag(C))
	}
	ni <- colSums(abs(key))
	#R <- cov2cor(C)
	r <- t(r)
	re <- r[,criteria]   
}
short.key <- list()
value <- list()
#R is the correlation with the criterion
#re is the correlations of each item with the criteria

for(i in 1:ny) {short.key[[criteria[i]]] <- round(short(key[,i,drop=FALSE],r),digits=digits) 

if(!is.null(dictionary)) {if(!is.factor(dictionary)) {temp <- lookup(rownames(short.key[[criteria[i]]]),dictionary)

  value[[criteria[[i]]]] <- merge(short.key[[i]],temp,by="row.names",all.x=TRUE,sort=FALSE)
  rownames( value[[criteria[[i]]]]) <-  value[[criteria[[i]]]][,1]
  value[[criteria[[i]]]] <- value[[criteria[[i]]]][-1]             #this looks weird but is because there is an extra name
  ord <- order(abs(value[[criteria[[i]]]][[criteria[[i]]]]),decreasing=TRUE)

  value[[criteria[[i]]]] <- value[[criteria[[i]]]][ord,]
 } 
 }}


results <- list(r=R,n.items=ni,R=re,cut=cut,short.key=short.key,value=value,key=key,ordered=ord.name,scores=score)
class(results) <- c("psych","bestScales")
return(results)
}

 
 
 "bestReplicatedItems" <- function( L) {
    n.iters <- length(L)
    n.vars <- length(L[[1]] )
    vars <- names(L[[1]])
    item.nums <- list()
    one.criterion <- list()
    for (j in 1:n.vars) { 
    for (i in 1:n.iters) {one.criterion[[i]] <- L[[i]][j] } 
    select <- sub("-","",unlist(one.criterion))
     item.nums[[vars[j]]] <- sort(table(select),decreasing=TRUE)
     }
     return(item.nums)
         }     
    
 
 print.psych.bestScales <- function(x,digits=2) {
    if(!is.null(x$first.result)) {
   cat("\nCall = ")
   print(x$Call)
   # print(x$first.result)
   #  print(round(x$means,2))
    print(x$summary,digits=digits)
      x$replicated.items
      
    items <- x$items
     size <- NCOL(items[[1]])
     nvar <- length(items)
     for(i in 1:nvar) {
     if(length(items[[1]]) >3 ) items[[i]] <- items[[i]][-1]
     items[[i]][1:3] <- round(items[[i]][1:3],digits)
      }

     cat("\n Best items on each scale with counts of replications\n")
     print(items)} else {
     df <- data.frame(correlation=x$r,n.items = x$n.items)
    cat("The items most correlated with the criteria yield r's of \n")
    print(round(df,digits=digits)) 
    if(length(x$value) > 0) {cat("\nThe best items, their correlations and content  are \n")
     print(x$value) } else {cat("\nThe best items and their correlations are \n")
     for(i in 1:length(x$short.key)) {print(round(x$short.key[[i]],digits=digits))} 
     }  
     } 
  }