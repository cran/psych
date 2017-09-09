#added  just correlate with criteria to speed it up (June 23, 2017)
"bestItems" <- 
function(x,criteria=1,cut=.3, abs=TRUE, dictionary=NULL,digits=2) {
if(!isCorrelation(x))  x <- cor(x,x[criteria],use="pairwise") #convert to correlation if necessary
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
 function(x,criteria,cut=.1,n.item =10,overlap=FALSE,dictionary=NULL,impute="none", n.iter =1,frac=.9,digits=2) {

  #first, define function to be parallelized
   #mcmapply for parallel, mapply for debugging 
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
	 if(n.iter > 1) {
  		cross <- scoreFast(key.list,x[-ss,],impute=impute,min=1,max=6)
 	 validity <- diag(cor(cross,x[-ss,criteria],use="pairwise"))
  	result <- list(validity=c(scores$r,validity),keys=key.list)
  } else {result <- scores
          result$key.list <- key.list}
 return(result) }
 
 
 n.obs <- nrow(x)
 if(n.iter ==1) {frac <- 1   #don't bother to parallelize, just do one trial
   result <- short(1,x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=frac)
   } else {  #the case for n.iter > 1.  We want to parallelize this because we are working pretty hard
 
result <- list()

result <- mcmapply(short,c(1:n.iter),MoreArgs=list(x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=frac))
  
 #save the keys and the summary 
  validity <- list()
  keys <- list()
  for(i in (1:n.iter)) { validity[[i]] <- result[[i*2 -1]]
    keys[[i]] <- result[[i*2]]}
  
   result.df <- data.frame(matrix(unlist(validity),ncol=2*length(criteria),byrow=TRUE))
  colnames(result.df) <-c(paste("derivation",criteria),paste("validation",criteria)) 
  result <- list(validity = result.df,keys = keys)}
  return(result)
}

"bScales" <- 
function(x,criteria,cut=.1,n.item =10, overlap=FALSE,dictionary=NULL,impute="median",digits=2) {

 #created 20/2/14
 #find the scales based upon the items that most correlate with a criteria
 #pure dust bowl empiricism
 #modified 13/3/15 to handle the problem of missing item labels
 #Completely redone June 2017 to allow for raw data 
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
 if(ny > 1 ) {ord <- apply(abs(r[,criteria]),2,order,decreasing=TRUE) 
     for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[ord[n.item[i]+1,i],criteria[i]])) 
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
 if(!overlap)  {key[criteria,criteria] <- 0} else {for(i in 1:ny) key[criteria[i],criteria[i]] <- 0}
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
  value[[criteria[[i]]]] <- value[[criteria[[i]]]][-1]
  ord <- order(abs(value[[criteria[[i]]]][[criteria[[i]]]]),decreasing=TRUE)

  value[[criteria[[i]]]] <- value[[criteria[[i]]]][ord,]
 } 
 }}

results <- list(r=R,n.items=ni,R=re,cut=cut,short.key=short.key,value=value,key=key,ordered=ord.name,scores=score)
class(results) <- c("psych","bestScales")
return(results)
}
 