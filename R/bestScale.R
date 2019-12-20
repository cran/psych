
#Modified October, 2019 to have a variable number of items per scale
#Modified Sept 8, 2019 to include weighted scoring  
#Modified April 15, 2019 to make structure somewhat cleaner


"bestScales" <- 
 function(x,  #the raw data or a correlation matrix
          criteria, # the criteria (name) to predict
          min.item=NULL,max.item=NULL,delta=0,  #parameters for multiple solutions 
          cut=.1,n.item =10,wtd.cut = 0,wtd.n=10, #just one solution / criteria
          n.iter =1,folds=1,p.keyed=.9, #how many bootstraps (n.iter) or folds
          overlap=FALSE,dictionary=NULL,check=TRUE,impute="none", log.p=FALSE, digits=2) {
 cl <- match.call()
  first <- TRUE
 #check for bad input   -- the Mollycoddle option 
	if(is.vector(criteria) & any( !(criteria %in% colnames(x)) )) {
 	cat("\nCriteria names are incorrectly specified. Offending items are ", criteria[which(!(criteria %in% colnames(x)))],"\n")
 	stop("I am sorry.  I am stopping because you did not specify the criteria correctly.  See above. ")} 
  #further check if the dictionary is specified correctly
  if(!is.null(dictionary)) if(length(dictionary) < 1) stop("I am sorry.  I am stopping because you did not specify the dictionary correctly.   ")
#check and delete variables with no variance (by default)
if(check) {item.var <- apply(x,2,sd,na.rm=TRUE)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0) ) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
             }
             }
  
#check various parameters
#
frac <- 1
if(folds > 1) {frac = 1/folds
if(n.iter !=folds) {n.iter <- folds
 cat('Number of iterations set to the number of folds = ',n.iter) }
 }
 set.seed(NULL)
 old.seed <- .Random.seed[42]   #we save this if we want to do k-fold cross validation
 
####
#first, define function to be parallelized
#mcmapply for parallel, mapply for debugging 

short <- function(i,x,n.obs,criteria,cut,n.item,impute,digits,dictionary,frac,log.p=FALSE,min.item,max.item) {#this is the function that is iterated 
		   multi.score <- NULL
           multi.cross <- NULL

	if(n.iter > 1) {
	   if(!isCorrelation(x)) { ss <- (1:n.obs) 
	  if(frac==1) {ss <- sample(ss,n.obs,replace=TRUE)  #bootstrap resampling  ss is 'in the bag'
	  }  else {
	  	set.seed(old.seed) #this will take the same random sequence for this run so that we can do k fold
		ss <- sample(ss,n.obs,FALSE)  # this is the 1:n.obs in the same random order for each of the k fold
		ss <- ss[-(((i-1)*frac*n.obs +1):(i*frac*n.obs))]     #this drops frac cases out each trial  
	    }
	   #the main function for finding the best items is right here
	   #probably don't need to pass dictionary every time, since we rarely examine a single pass
	 	scores   <-  bScales(x[ss,],criteria=criteria,cut=cut,
	             n.item =n.item,overlap=overlap,dictionary=dictionary,impute=impute,digits=digits,log.p=log.p) 
	             
	          #These next two lines then try to find the optimal number of items for this pass  
	          #Not clear if we really need to do this for every iteration, perhaps just for the final pooled r
	            if(!is.null(min.item)){

	             multi.score <- fastValidity(items=x[ss,],criteria=criteria,r=NULL,overlap=overlap, nlow=min.item,nhigh=max.item)
	             multi.cross <- fastCrossValidity(new.items = x[-ss,],r=multi.score$r,item.order=multi.score$item.order,criteria=criteria,
	                nlow=min.item,nhigh=max.item,delta=0,overlap=overlap,optimal.n=multi.score$optimal.unit.n,optimal.wtd.n = multi.score$optimal.wtd.n,)} else {multi.score <- NULL} 
	             } else {message("iterative solutions not possible for correlation matrices")
	            n.iter <- 1 
	            }} else { # a correlation matrix or n.iter = 1
		scores   <-  bScales(x,criteria=criteria,cut=cut,
	             n.item =n.item,overlap=overlap,dictionary=dictionary,impute=impute,digits=digits,log.p=log.p)
	             if(!is.null(min.item)){ multi.score <- fastValidity(items=x,criteria=criteria,r=NULL,overlap=overlap, nlow=min.item,nhigh=max.item)} else {multi.score <- NULL} 
	             }
      
	 key.list <- keys2list(scores$key) #this converts the -1 and 1s to a list with the variable names
	  	
	 if(n.iter > 1) {
  		cross <- scoreFast(key.list,x[-ss,],impute=impute,min=1,max=6) #why are these fixed values?  
 	    validity <- diag(cor(cross,x[-ss,criteria],use="pairwise"))
 	 #now, add the two new functions FastValidity and FastCrossValidity  
 	 
   #if we just want to do the optimal number of items on the summaries, we don't nee to return the multi.scores here
  	short.result <- list(r = c(scores$r,validity),key.list=key.list,R = scores$R,multi.scores=multi.score,multi.cross=multi.cross)
  } else {short.result <- scores    # this is the list of various objects from bScales
          short.result$key.list <- key.list
          short.result$multi.score <- multi.score
          short.result$multi.cross <- multi.cross}
     class(short.result) <- cbind("psych","bestScales")
     
 return(short.result) }  #this is the result from 1 iteration of all criteria
 ###
 ###
 
#begin the main function 
#if criteria is a separate data frame, combine x and criteria
#there are several alternative forms for criteria
#it is either a column name of x, or it is a separate data.frame/matrix
if(!is.null(dim(criteria))| (length(criteria) == NROW(x)))  { x <- cbind(x,criteria)
    if(NCOL(criteria) > 1 ){criteria <- colnames(criteria) }  else {criteria <- "criteria"}
  #criteria <- colnames(criteria)
    }
 
 n.obs <- nrow(x)
 #if((n.iter ==1)| first ) {   #don't bother to parallelize, just do one trial
 if((n.iter ==1)) { 
   first.result <- short(1,x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=1,min.item=min.item,max.item=max.item)
   first <- FALSE
   result <- first.result
   } else {first.result <- NULL}
     #the case for n.iter > 1.  We want to parallelize this because we are working pretty hard

 if(n.iter > 1) { 
result <- list()
#This does the work across n.iter and across all criteria
result <- mcmapply(short,c(1:n.iter),MoreArgs=list(x,n.obs=n.obs,criteria=criteria,cut=cut,n.item=n.item,impute=impute,digits=digits,dictionary=dictionary,frac=frac,min.item=min.item,max.item=max.item))

#we have done the heavy lifting, now we need to prepare various results for output.
if(delta >  0) { delta <- delta /sqrt(n.obs)}
result <- organize.results(result,x,n.iter=n.iter,p.keyed=p.keyed,dictionary=dictionary,wtd.cut=wtd.cut,wtd.n = wtd.n,overlap=overlap,min.item=min.item,max.item=max.item,delta=delta)  #makes the function a bit cleaner by doing this in its own function
#save the keys and the summary 

   
  } else {  #we just did a single pass, the simple summaries are already there
   result$best.keys=result$key.list

   final.means <-  colMeans(result$scores,na.rm=TRUE)
   final.sd <- apply(result$scores,2,sd,na.rm=TRUE)
 	if(length(criteria) > 1 ) {crit.mean <- colMeans(x[,criteria],na.rm=TRUE)
 	crit.sd <- apply(x[,criteria],2,sd,na.rm=TRUE)} else {
 	crit.mean <- mean(x[,criteria],na.rm=TRUE)
 	crit.sd <- sd(x[,criteria],na.rm=TRUE)}
   	result$final.stats <- data.frame(mean=final.means,sd=final.sd,r=result$r,crit.m=crit.mean,crit.sd =crit.sd)
    result$items <- NULL
    
   }
    result$Call <- cl
   result$first.result  <- first.result
  class(result) <- c("psych","bestScales")
  return(result)
}
#####################


#######################################
#This function takes the results from short for many trials and then tries to make sense of them
######################################
organize.results <- function(result,x=NA,n.iter=1,p.keyed=.9,dictionary=NULL,wtd.cut,wtd.n,overlap=overlap, min.item=min.item,max.item=max.item,delta=delta) {  
#The results are n.iter lists, each with validity,keys,R, and the results from multi.score

  validity <- list()
  #validity is a list of  elements repeated n.iter times
  #first are the validities
  #then are the keys
  #then are the item by criteria correlations
  #then the multi.score matrices
  keys <- R.list <- multi.valid <- multi.cross <-  list()
 
  
 #take the list from all the iterations, and organize them in a more meaningful way
  for(i in (1:n.iter)) { 
      validity[[i]] <- result[["r",i]]
       keys[[i]] <- result[["key.list",i]]
      R.list[[i]] <- result [["R",i,drop=FALSE]]
 if(!is.null(min.item)) {     multi.valid [[i]] <- result[["multi.scores",i,drop=FALSE]] 
       multi.cross [[i]] <- result[["multi.cross",i,drop=FALSE]] 
       }
    }
     
    replicated.items <- bestReplicatedItems(keys)   

   items <- list()
    item.mean <- list()
    best.keys <- list()
   criteria <- names(replicated.items)
   optimal.n <- optimal.wtd.n <- optimal.unit.deriv <- optimal.wtd.deriv <- optimal.cross.unit <- optimal.cross.wtd <- cross.n <- cross.wtd.n <- list()
   
   #we can find the optimal length for all criteria at once
   if(!is.null(min.item)) {
   for (i in 1:n.iter) {
     optimal.n[[i]] <- multi.valid[[i]][["optimal.unit.n"]]
     optimal.wtd.n[[i]] <- multi.valid[[i]][["optimal.wtd.n"]]
     optimal.unit.deriv[[i]] <-  multi.valid[[i]][["optimal.unit.deriv"]]
     optimal.wtd.deriv[[i]] <- multi.valid[[i]][["optimal.wtd.deriv"]]
     optimal.cross.unit[[i]] <-  multi.cross[[i]][["cross.unit"]]
     optimal.cross.wtd[[i]] <-  multi.cross[[i]][["cross.wtd"]]
     cross.n[[i]] <-  multi.cross[[i]][["cross.n"]]
     cross.wtd.n[[i]] <-  multi.cross[[i]][["cross.wtd.n"]]
     }
     optimal.n.mean <-apply(matrix(unlist(optimal.n),nrow=n.iter,byrow=TRUE),2,median) 
     optimal.wtd.mean <- apply(matrix(unlist(optimal.wtd.n),nrow=n.iter,byrow=TRUE),2,median) 
     optimal.unit.deriv <- colMeans(matrix(unlist(optimal.unit.deriv),nrow=n.iter,byrow=TRUE)) 
     optimal.wtd.deriv <- colMeans(matrix(unlist(optimal.wtd.deriv),nrow=n.iter,byrow=TRUE)) 
     optimal.cross.unit <- colMeans(matrix(unlist(optimal.cross.unit),nrow=n.iter,byrow=TRUE)) 
     optimal.cross.wtd <- colMeans(matrix(unlist(optimal.cross.wtd),nrow=n.iter,byrow=TRUE)) 
     cross.n <- apply(matrix(unlist(cross.n),nrow=n.iter,byrow=TRUE),2,median) 
     cross.wtd.n <- apply(matrix(unlist(cross.wtd.n),nrow=n.iter,byrow=TRUE),2,median) 

}

  
   #but we need to find item statistics one criteria at a time
   for(j in 1:length(criteria)) {
     #first, find  the means and standard deviations for each selected item
     rep.item <-  replicated.items[[j]][replicated.items[[j]] >= n.iter * p.keyed]
     if(length(rep.item)==0) rep.item <- replicated.items[[j]][1]
   #   if(length(criteria) > 1 ) {for (i in 1:n.iter) { item.mean[[i]] <-  R.list[[i]][names(replicated.items[[j]][replicated.items[[j]] > n.iter * p.keyed]),criteria[j]] }
    #  } else { for (i in 1:n.iter) {item.mean[[i]] <- R.list[[i]][names(replicated.items[[j]][replicated.items[[j]] > n.iter * p.keyed])] } }
  
  
  
   for (i in 1:n.iter) {if(length(criteria) > 1) {
   item.mean[[i]] <-  R.list[[i]][names(rep.item),criteria[j]]} else {item.mean[[i]] <-  R.list[[i]][names(rep.item)]}
   }

     item.m <- matrix(unlist(item.mean),nrow=n.iter,byrow=TRUE)
     colnames(item.m) <- names(rep.item)
     means = colMeans(item.m,na.rm=TRUE)
     sds <- apply(item.m,2,sd,na.rm=TRUE)  
  #   Freq <- colSums(!is.na(item.m))  #This is the total number of items and just reflect n.iter
      Freq <- as.vector(rep.item)
      names(Freq) <- names(rep.item)

#items [[criteria[j] ]] <-  cbind(replicated.items[[j]],Freq=Freq,mean.r=means,sd.r = sds,dictionary[names(replicated.items[[j]]),])
    items[[criteria[j]]] <-  cbind(Freq,mean.r = means,sd.r = sds,dictionary[names(rep.item),])
	items[[criteria[j]]] <- psychTools::dfOrder(items [[criteria[j] ]],"-mean.r",absolute=TRUE) #sort on the mean.r column
#	items[[criteria[j]]] <- items[[criteria[j]]][items[[criteria[j]]][,"Freq"] >= n.iter * p.keyed,]
	
	#now prepare the best.keys list

   if(!is.null(dim(items[[criteria[[j]] ]] ))){ direction <- sign(items[[criteria[[j]] ]][,"mean.r"]) 
         direction <- as.matrix(direction)
         rownames(direction) <- rownames(items[[criteria[[j]] ]])
         count <- items[[criteria[[j]]]][,"Freq"]} else {
        if(!is.null(items[[criteria[[j]] ]])) {
       # items [[criteria[j] ]] <-  cbind(Freq=Freq,mean.r=means,sd.r = sds,dictionary[names(replicated.items[[j]]),])
       items [[criteria[j] ]] <-  cbind(Freq=Freq,mean.r=means,sd.r = sds,dictionary[names(replicated.items[[j]]),])
         direction <-  sign(items[[criteria[[j]] ]]["mean.r"]) 
              names(direction) <- names(Freq)
               direction <- as.matrix(direction)
        
              count <- items[[criteria[[j]]]][1] }
                else {count <- 0}
           }
        count <- count >= n.iter*p.keyed   
         
   
   if(sum(count,na.rm=TRUE) > 0) {
    best.keys[[j]] <- rownames(direction)[count]
    direction <- direction[count,drop=FALSE]
     if(length(direction)> 1) { best.keys[[j]][direction < 0] <- paste0("-", best.keys[[j]][direction < 0]) }
     if((length(direction) ==1) && (!is.na(direction))) {
   best.keys[[j]][direction < 0] <- paste0("-", best.keys[[j]][direction < 0]) }
   } else { best.keys[[j]] <- NA
}
    }
    
#Find the mean, zero order correlation of each item with each criteria
#We do this by pooling the data in R.list
mean.raw.r <- matrix(unlist(R.list),ncol=NCOL(R.list[[1]]) * NROW(R.list[[1]]),byrow=TRUE )
sd.raw.r <- apply(mean.raw.r,2,sd,na.rm=TRUE)
sd.raw.r <- matrix(sd.raw.r,ncol=length(criteria))
mean.raw.r <- matrix(colMeans(mean.raw.r,na.rm=TRUE),ncol=length(criteria))
if(length(criteria) == 1) {colnames(mean.raw.r) <- criteria
  rownames(mean.raw.r) <- names(R.list[[1]])}  else {colnames(mean.raw.r) <- colnames(R.list[[1]])
rownames(mean.raw.r) <- rownames(R.list[[1]])}
final.mean.r <- mean.raw.r
mean.raw.r[abs(mean.raw.r) < wtd.cut] <- 0
#now, drop all except the wtd.n items


 
ny <- length(criteria)
nvar <- NROW(mean.raw.r)
if(ny > 1 ) {ord <- apply(abs(mean.raw.r[,criteria]),2,order,decreasing=TRUE) 
     for (i in 1:ny) {mean.raw.r[ord[(wtd.n+1):nvar,i],criteria[i]]  <- 0 
    }
     } else {
         ord <- order(abs(mean.raw.r),decreasing=TRUE)
         for (i in 1:ny) {mean.raw.r[ord[(wtd.n+1):nvar]] <- 0 }
         }
        
N.wtd <- colSums(abs(mean.raw.r) >0)
    
   if(length(best.keys) == length(criteria)) names(best.keys) <- criteria 	

   #Find the results for best keys	
   	final.scale <- scoreFast(best.keys,x)   #these are the unit weighted
   	final.raw.scale <- scoreWtd(mean.raw.r,x)   #these are the zero order weighted scores
   	final.raw.valid <- diag(cor(final.raw.scale,x[,criteria,drop=FALSE],use="pairwise") )
   final.valid <- diag(cor(final.scale, x[,criteria,drop=FALSE],use="pairwise")	)
   final.means <-  colMeans(final.scale,na.rm=TRUE)
   final.sd <- apply(final.scale,2,sd,na.rm=TRUE)
   crit.mean <- colMeans(x[,criteria,drop=FALSE],na.rm=TRUE)
   crit.sd <- apply(x[,criteria,drop=FALSE],2,sd,na.rm=TRUE)
   
   
   
 
   result.df <- data.frame(matrix(unlist(validity),ncol=2*length(criteria),byrow=TRUE))
  colnames(result.df) <-c(paste("derivation",criteria),paste("validation",criteria)) 

  if(!is.null(min.item)) {
  multi.derivation.df <- data.frame(n=optimal.n.mean,unit=optimal.unit.deriv,n.wtd=optimal.wtd.mean,wtd=optimal.wtd.deriv,valid.n=cross.n,valid.unit=optimal.cross.unit,valid.wtd.n = cross.wtd.n,valid.wtd=optimal.cross.wtd)
  rownames(multi.derivation.df ) <- criteria
   } else {multi.derivation.df <- NULL}
ncriteria <- length(criteria)
 if(!is.null(min.item)){ 

 	final.multi.validities <- fastValidity(x,criteria,r=final.mean.r, nlow=min.item, nhigh=max.item,overlap=overlap)
 	final.order <- final.multi.validities$item.order
 	final.item.valid.list <- list()
  	for (j in 1 : ny ) {if(!is.null(dictionary)) {final.item.valid.list[[criteria[j]]] <- data.frame(item=rownames(final.mean.r)[final.order[1:max.item,j]]
  	        ,r=final.mean.r[final.order[1:max.item,j],j],unit = final.multi.validities$unit.deriv[,j],wtd=final.multi.validities$wtd.deriv[,j],dictionary[rownames(final.mean.r)[final.order[1:max.item,j]],])
    	 } else { final.item.valid.list[[criteria[j]]] <- data.frame(item=rownames(final.mean.r)[final.order[1:max.item,j]],r=final.mean.r[final.order[1:max.item,j],j],unit = final.multi.validities$unit.deriv[,j],wtd=final.multi.validities$wtd.deriv[,j])}
  }
  
  } else {final.multi.validities <- NULL
    final.item.valid.list <- NULL}
  
#now, organize the output object into a reasonable order
result <- list()
results <- list()  #to hold things we don't actually want to return

 #now get out the items and incremental validities for each scale
  results$means = colMeans(result.df,na.rm=TRUE)
  results$sd <- apply(result.df,2,sd,na.rm=TRUE)
  
result$summary <- data.frame(derivation.mean= results$means[1:ncriteria],derivation.sd = results$sd[1:ncriteria],validation.m=results$mean[(ncriteria+1):(2*ncriteria)],
            validation.sd =results$sd[(ncriteria+1):(2*ncriteria)],final.valid = final.valid,final.wtd=final.raw.valid,N.wtd=N.wtd )
 rownames(result$summary) <- criteria   
         
#result <- list(validity = result.df,multi.validities=multi.derivation.df,items=items,replicated.items =replicated.items,keys = keys,final.mean.r =final.mean.r,multi.validities=multi.valid)
 
 
 result$optimal <- multi.derivation.df
 result$best.keys <- best.keys
 result$weights <- mean.raw.r
 result$final.item.list <- final.item.valid.list
 result$multi.validities <- final.multi.validities
 result$items <- items
 if(!is.null(min.item)) {
    result$optimal.keys <- optimal.keys(final.multi.validities,delta=delta)
    result$optimal.weights <- optimal.weights(final.multi.validities,delta=delta)
 n.optimal.unit <- sapply(result$optimal.keys,length)
 n.optimal.wtd <-  apply(result$optimal.weights,2,function(x) sum(abs(x) > 0) )
 result$optimal <- data.frame(result$optimal,n.final=n.optimal.unit,n.wtd.final = n.optimal.wtd)
 }
  result$stats <- data.frame(mean=results$means,se=results$sd)
  
   # result$final <- final.valid
    result$scores <- final.scale
    result$wtd.scores <- final.raw.scale
    result$final.stats <- data.frame(mean=final.means,sd=final.sd,r=final.valid,crit.mean = crit.mean,crit.sd=crit.sd,final.wtd=final.raw.valid,N.wtd=N.wtd)
   
   # result$sd.weights <- sd.raw.r
   # result$final.raw <- final.raw.valid
    return(result)
   } 
 ###########################################  #end of organize results
 ########################################### 
 
 
  #This one actually does the work  -- but should not process n.item at this point
"bScales" <- 
function(x,criteria,cut=.1,n.item =10, overlap=FALSE,dictionary=NULL,impute="median",digits=2,log.p=FALSE) {

 #created 20/2/14
 #find the scales based upon the items that most correlate with a criteria
 #pure dust bowl empiricism
 #modified 13/3/15 to handle the problem of missing item labels
 #Completely redone June 2017 to allow for raw data and bootstrapping
 ##



#begin the main function
 ##Basically two cases:
 #a correlation matrix is provided and we do basic matrix algebra
 #or raw data are provided (getting ready to do bootstrapping) and we find just the necessary correlations
 
 nvar <- ncol(x)
 if(isCorrelation(x)) {r <- x      #  case 1
    raw <- FALSE} else {  #case 2
    y <- x[,criteria]
    if(log.p) {r <- log(corr.test(x,y)$p)} else { r <- cor(x,y,use="pairwise")}
   
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
if(length(n.item) == 1) n.item <- rep(n.item,ny) #

#We have the correlations with the criteria, we can 

#this next part just finds the cut values to use
 if(!overlap)  {r[criteria,criteria] <- 0} else {for(i in 1:ny) r[criteria[i],criteria[i]] <- 0}
 if(ny > 1 ) {ord <- apply(abs(r[,criteria]),2,order,decreasing=TRUE) 
     for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[ord[n.item[i],i],criteria[i]])) 
    }
     } else {
         ord <- order(abs(r[,criteria]),decreasing=TRUE)
         for (i in 1:ny) {cut[i] <- max(cut[i],abs(r[ord[n.item[i]+1],criteria])) }
        }
#    cut has been adjusted

#The unit weights
 key <- matrix(0,ncol=ny,nrow=nvar)
 key[t(t(r[,criteria]) >= cut)] <- 1
 key[t(t(r[,criteria]) <= -cut)]<- -1
 rownames(key)  <- rownames(r)
 colnames(key)  <- criteria

 k <- key  #this just gets it to be a matrix of the right size and names

 #colnames(key) <- paste(criteria,"S",sep=".")
 #colnames(key) <- criteria
 #now, drop those items from the keys that are not used
 used <- rowSums(abs(key))
 key <- key[used > 0,,drop=FALSE]  
 x <- x[,used >0,drop=FALSE]
 
#now, if we have raw data, find the correlation of the composite scale with the criteria
#if we have raw data, then we find the scales from the data 
if(raw)  { #case 2
#score <- matrix(NA,ncol=ny,nrow=nrow(x))
#for (i in (1:ny)) {
#   score[,i] <- rowSums(t((key[,i])* t(x)),na.rm=TRUE)
 #     }
 score <- scoreFast(key,x,impute=impute,min=1,max=6)     #min and max should not be fixed values
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

for(i in 1:ny) {short.key[[criteria[i]]] <- round(key.value(key[,i,drop=FALSE],r),digits=digits) 
 #actually we should not bother with the dictionary here, just at the summary level 
if(!is.null(dictionary)) {if(!is.factor(dictionary)) {temp <- lookup(rownames(short.key[[criteria[i]]]),dictionary)

  value[[criteria[[i]]]] <- merge(short.key[[i]],temp,by="row.names",all.x=TRUE,sort=FALSE)
  rownames( value[[criteria[[i]]]]) <-  value[[criteria[[i]]]][,1]
  value[[criteria[[i]]]] <- value[[criteria[[i]]]][-1]             #this looks weird but is because there is an extra name
  ord <- order(abs(value[[criteria[[i]]]][[criteria[[i]]]]),decreasing=TRUE)

  value[[criteria[[i]]]] <- value[[criteria[[i]]]][ord,]
 } 
 }}

bScales.results <- list(r=R,n.items=ni,R=re,cut=cut,short.key=short.key,value=value,key=key,ordered=ord.name,scores=score)
class(bScales.results) <- c("psych","bestScales")      #This is the solution for one pass
return(bScales.results)
}

################################
 
 #various minor functions used in bestScales  

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

key.value <- function(key,r) { 
	 kn <- names(key[abs(key[,1]) > 0,1])
	if(is.null(kn)) kn <- names(which(abs(key[,1]) > 0))
	 cn <- colnames(key)
 	ord <- 	order(abs(r[kn,cn]),decreasing=TRUE)
 	kn <- kn[ord]
	 result <- r[kn,cn,drop=FALSE]
	 return(result)
	}

#


"finalsummary" <- function(r,keys) {
  }
 
 
 "bestReplicatedItems" <- function( L) {
    n.iters <- length(L)
   # n.vars <- NCOL(L[[1]]) 
    vars <- names(L[[1]])
    n.vars<- length(vars)
    item.nums <- list()
    one.criterion <- list()
    for (j in 1:n.vars) {  #do this over the criteria
    for (i in 1:n.iters) {one.criterion[[i]] <- L[[i]][j] } 
    select <- sub("-","",unlist(one.criterion))
     item.nums[[vars[j]]] <- sort(table(select),decreasing=TRUE)
     nam.items <- names(item.nums)
     }
     item.nums <- as.vector(item.nums)
     names(item.nums) <- nam.items
     return(item.nums)
         }     
    
predict.bestScales <- function(object,data,new.data) 
{keys <- object$keys
if (is.null(keys)){ keys<- object$key.list
 keys <- make.keys(data,keys) }  
stats <- describe(data,fast=TRUE,ranges=FALSE)
old.mean <- stats$mean
old.sigm <- stats$sd
z.data <- scale(new.data,center=stats$mean,scale=stats$sd)
z.data[is.na(z.data)] <- 0
predicted <-  z.data   %*% keys 
return(predicted)
}


#does not do what I want
predict.wtdScales <- function(object,data,new.data) 
{weights <- object$weights
  
stats <- describe(data,fast=TRUE,ranges=FALSE)
old.mean <- stats$mean
old.sigm <- stats$sd
predicted <- scoreWtd(weights,new.data)
predicted <- scale(predicted,center=stats$mean,scale=stats$sd)
return(predicted)
}

#####
#print.psych.bestscales is called from the psych.print function

print.psych.bestScales <- function(x,digits=2,short=NULL,...) {
   cat("\nCall = ")
   print(x$Call)
if(!is.null(x$items)) {
    print(x$summary,digits=digits)
    if(!is.null(x$optimal)) {
     cat("\n Optimal number of items, derivation and cross validation\n")
   print(x$optimal,digits=digits) }
     # x$replicated.items
      
    items <- x$items
     size <- NCOL(items[[1]])
     nvar <- length(items)
 if(is.null(x$optimal)) {
   cat("\n Best items on each scale with counts of replications\n")    
 for(i in 1:nvar) {
 cat("\n Criterion = ",names(items[i]),"\n")
  if(length(items[[i]]) > 3) {
 temp <- data.frame(items[[i]])
   temp[2:3] <- round(temp[2:3],digits)} else{temp <- items[[i]]
   temp[2:3] <- round(temp[2:3],digits)
  }
   print(temp)  
 }} else {
 items <- x$final.item.list
  cat("\n Best items on each scale with cumulative validities\n")    
 for(i in 1:nvar) {
 cat("\n Criterion = ",names(items[i]),"\n")
  temp <- data.frame(items[[i]])
  temp <- temp[-1]
  if(length(items[[i]]) > 3) {
   temp[1:3] <- round(temp[1:3],digits)} else{temp <- items[[i]]
   temp[1:3] <- round(temp[1:3],digits)
  }
   print(temp)  
 }}
  
    # print(items)
     } else {
     df <- data.frame(correlation=x$r,n.items = x$n.items)
    cat("The items most correlated with the criteria yield r's of \n")
    print(round(df,digits=digits)) 
    if(length(x$value) > 0) {cat("\nThe best items, their correlations and content  are \n")
     print(x$value) } else {cat("\nThe best items and their correlations are \n")
     for(i in 1:length(x$short.key)) {print(round(x$short.key[[i]],digits=digits))} 
     }  
     } 
 }
#end of print.psych.bestScales        
      

    
#These next two functions were added in October, 2019 to allow for  finding the maximum value of cross.validated as a function of n.item  
#ValidityList <- mapply(FastValidity,c(1:ny),MoreArgs=list(nlow=nlow,nhigh=nhigh))  #probably not that helpful

### Find predictions and validities for scales from nlow to nhigh number of items 
fastValidity <- function(items, criteria,r=NULL,nlow,nhigh,overlap) { #r and criteria are found before
 ny <-length(criteria)
 if(is.null(r)) r <- cor(items,items[criteria],use="pairwise")
wtd.validity <-  unit.validity <- matrix(NA,nrow=nhigh,ncol= ny)
if(!overlap)  {r[criteria,criteria] <- 0} else {for(i in 1:ny) r[criteria[i],criteria[i]] <- 0} 

colnames(unit.validity) <- criteria 
colnames(wtd.validity) <- criteria

item.min <- min(items,na.rm=TRUE)
item.max <- max(items,na.rm = TRUE)
if(item.max < 10) {
item.range.correction <- item.max - item.min + 1} else item.range.correction <- 7   #this is the case where we include some weird items, like age
for(scale in 1:ny) {
 if(ny > 1 ) {ord <- apply(abs(r[,criteria,drop=FALSE]),2,order,decreasing=TRUE) 
     } else {
         ord <- matrix(order(abs(r[,criteria,drop=FALSE]),decreasing=TRUE))
              rownames(ord) <- rownames(r)
              colnames(ord) <- criteria
        }
    
         abs.item <- t(t(items) * sign(r[,scale]) + sign(r[,scale] < 0) *item.range.correction )   #this gets all the items to be scored and adds in the max - min + 1       
         wtd.item <- t(t(items) * r[,scale] + sign(r[,scale] < 0) *item.range.correction  )    #these are the wtd item scores 
    for (j in nlow: nhigh){
    #  temp <- abs.item[,ord[1:j,scale,drop=FALSE]] 
    #  wtd.temp <- wtd.item[,ord[1:j,scale,drop=FALSE]]      	
    if(j > 1) {scores <- rowMeans( abs.item[,ord[1:j,scale,drop=FALSE]] ,na.rm=TRUE)
              wtd.scores <- rowMeans(wtd.item[,ord[1:j,scale,drop=FALSE]],na.rm=TRUE)
             } else {scores <- abs.item[,ord[1:j,scale,drop=FALSE]]
                     wtd.scores <- wtd.item[,ord[1:j,scale,drop=FALSE]]}
    unit.validity[j,scale] <- cor(scores,items[,criteria[scale]],use="pairwise")
    wtd.validity [j,scale] <- cor(wtd.scores,items[,criteria[scale]],use="pairwise") 
    }
  }
  optimal.unit.n <- apply(unit.validity,2,which.max)
  optimal.wtd.n <- apply(wtd.validity,2,which.max)
  optimal.unit.valid <- apply(unit.validity,2,max)
  optimal.wtd.valid <- apply(wtd.validity,2,max)
  result <- list(optimal.unit.n=optimal.unit.n,  optimal.wtd.n = optimal.wtd.n,  
                 optimal.unit.deriv=optimal.unit.valid, optimal.wtd.deriv=optimal.wtd.valid,
                 unit.deriv=unit.validity,wtd.deriv = wtd.validity,item.order = ord,r=r,item.range.correction=item.range.correction)
 return(result) 
  } 
  
 
 #This takes the item order from fastValidity and applies it to a new data set,using the old correlations 
fastCrossValidity <- function(new.items,r,item.order,criteria,nlow,nhigh,overlap,optimal.n,optimal.wtd.n,delta=0,item.range.correction=0) { #r  and order are from the derivation set
 ny <-length(criteria)
   #r and item.order are from the derivation sample
wtd.cross <-  unit.cross <- matrix(NA,nrow=nhigh,ncol= ny)
if(!overlap)  {r[criteria,criteria] <- 0} else {for(i in 1:ny) r[criteria[i],criteria[i]] <- 0} 

colnames(unit.cross) <- criteria 
colnames(wtd.cross) <- criteria
ord <- item.order

for(scale in 1:ny) {
 
         abs.item <- t(t(new.items) * sign(r[,scale]) + sign(r[,scale] < 0) * item.range.correction) 
   #this gets all the items to be scored        
         wtd.item <- t(t(new.items) * r[,scale] + item.range.correction)    #these are the wtd item scores 
    for (j in nlow: nhigh){ 
       temp <- abs.item[,ord[1:j,scale,drop=FALSE]] 
      wtd.temp <- wtd.item[,ord[1:j,scale,drop=FALSE]]       	
    if(j > 1) {scores <- rowMeans(temp[,1:j],na.rm=TRUE)
              wtd.scores <- rowMeans(wtd.temp[,1:j],na.rm=TRUE)
             } else {scores <- abs.item[,ord[1:j,scale,drop=FALSE]]
                     wtd.scores <- wtd.item[,ord[1:j,scale,drop=FALSE]]}
    unit.cross[j,scale] <- cor(scores,new.items[,criteria[scale]],use="pairwise")
    wtd.cross [j,scale] <- cor(wtd.scores,new.items[,criteria[scale]],use="pairwise") 
    }
   
 cross.unit.valid <- apply(unit.cross,2,max)
  cross.wtd.valid <- apply(wtd.cross,2,max)
  
 # temp <- apply(unit.cross,2,function(x) which(x >= (max(x) - delta)))
 #if(is.list(temp)) {cross.unit.n <- sapply(temp,function(x) x[1],simplify=TRUE) } else {cross.unit.n <- temp[1]} 
  cross.unit.n <- apply(unit.cross,2,which.max)
  cross.wtd.n <- apply(wtd.cross,2,which.max)
  
    
  optimal.cross.unit <- diag(unit.cross[optimal.n,1:ny,drop=FALSE])
  optimal.cross.wtd <-diag( wtd.cross[optimal.wtd.n,1:ny,drop=FALSE])
  }

  result <- list(unit.cross=unit.cross, wtd.cross = wtd.cross, cross.unit= optimal.cross.unit,
                  cross.wtd=optimal.cross.wtd, cross.n=cross.unit.n, cross.wtd.n=cross.wtd.n)
 return(result) 
  }  
  
optimal.keys <- function(L,delta=0) {
   #take the information from multi.validities and create a keys.list and a weights matrix
    criteria <- names(L[["optimal.unit.n"]])
    n <-L[["optimal.unit.n"]]
    unit.cross <- L[["unit.deriv"]]
    if(delta>0) {
    temp <- apply(unit.cross,2,function(x) which(x >= (max(x,na.rm=TRUE) - delta)))
 if(is.list(temp)) {n <- sapply(temp,function(x) x[1],simplify=TRUE) } else {if (is.vector(temp)) {n <- temp}  else {n <- temp[1,,drop=FALSE]}}
    } 
    item.order <- L [["item.order"]]
    var.names <- rownames(L[["r"]])
    r <- L [["r"]]
    keys <- direction <-  list()
    for (j in 1:length(criteria)) {   
      keys[[j]] <- var.names[item.order[1:n[j],j]]
      direction <- sign(r[item.order[1:n[j]],j])
      keys[[j]][direction <0 ] <- paste0("-",keys[[j]][direction < 0])
       }
       names(keys) <- criteria
       return(keys)
   }
   
   optimal.weights <- function(L,delta=0) {
   #take the information from multi.validities and create a keys.list and a weights matrix
    criteria <- names(L[["optimal.unit.n"]])
    n <-L[["optimal.unit.n"]]
    wtd.cross <- L[["wtd.deriv"]]
    if(delta>0) {
    temp <- apply(wtd.cross,2,function(x) which(x >= (max(x,na.rm=TRUE) - delta)))
 if(is.list(temp)) {n <- sapply(temp,function(x) x[1],simplify=TRUE) } else {if (is.vector(temp)) {n <- temp}  else {n <- temp[1,,drop=FALSE]}}
    } 
    item.order <- L [["item.order"]]
    var.names <- rownames(L[["r"]])
    weights <- L [["r"]]
    
    cut <- abs(diag(weights[diag(item.order[n,,drop=FALSE]),,drop=FALSE]))
    weights[t(abs(t(weights)) < cut)] <- 0
      return(weights)
   }
  