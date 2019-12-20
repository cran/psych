#created 7/12/16 
#just score items without a lot of stats
#basically scoreItems with all  the stats removed\
#Parallelized July 28, 2018 and report the number of responses/scale
#added the "mollycoddle" feature March 19, 2019 to help the clueless user
"scoreFast"  <-
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE, impute="none",delete=TRUE,  min=NULL,max=NULL,count.responses=FALSE,digits=2) {
 
 smallFunction <- function(scale,keys) {
      if(is.null(keys)) return(NULL)
          	pos.item <- items[,which(keys[,scale] > 0)]
            neg.item <- items[,which(keys[,scale] < 0)]
          	neg.item <- max + min - neg.item
           	sub.item <- cbind(pos.item,neg.item)
           	if(count.responses) rs <- rowSums(!is.na(sub.item))
   if(totals) { scores <- rowSums(sub.item,na.rm=TRUE)} else {scores <- rowMeans(sub.item,na.rm=TRUE) } 
    if(count.responses) {return(c(scores,rs))} else {return(scores) }
  }
  
   cl <- match.call()
   if(is.data.frame(keys)) stop("I think you reversed keys and items.  I am stopping")
   raw.data <- TRUE
   if(impute == FALSE)  impute <- "none"
   if(is.list(keys)) {select <- sub("-","",unlist(keys))
      select <- select[!duplicated(select)]
      select <- select[!is.na(select)]
      #check for bad input   -- the Mollycoddle option 
if(any( !(select %in% colnames(items)) )) {
 cat("\nVariable names are incorrect. Offending items are ", select[which(!(select %in% colnames(items)))],"\n")
 stop("Improper input.  See above. ")}
       } else {
      keys <- keys2list(keys)
        select <- selectFromKeyslist(colnames(items),keys)
      select <- select[!duplicated(select)]
      select <- select[!is.na(select)]}    #added 11/23/18
      items <- items[,select,drop=FALSE]
      keynames <- colnames(keys)
      keys <- make.keys(items,keys)   #added 9/9/16 
   keys <- as.matrix(keys)   #just in case they were not matrices to start with
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
     abskeys <- abs(keys)
     keynames <- colnames(keys)
   num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
  num.ob.item <- num.item   #will be adjusted in case of impute = FALSE
    if (!missing) items <-  na.omit(items) 
    n.subjects <- dim(items)[1]
    
   items <- as.matrix(items)
    
   # response.freq <- response.frequencies(items)
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
        
        }  else { #handle the case of missing data without imputation
           scores <- matrix(NaN,ncol=n.keys,nrow=n.subjects)

   scoresList <- mcmapply(smallFunction,c(1:n.keys),MoreArgs=list(keys=keys))  #the parallelized function
    
    }      
    
    if (is.null(ilabels)) {
    	if (totals) {#ilabels<- paste("S",1:n.keys,sep="")} else {
    	             #ilabels <- paste("A",1:n.keys,sep="")} }
  ilabels<- paste(keynames,"S",sep="-")} else {
    	             ilabels <- paste(keynames,"A",sep="-")} }
                    
 if(count.responses) { scores <- scoresList[1:n.subjects,]
   responses <- scoresList[(n.subjects+1):nrow(scoresList),]
   colnames(scores) <- ilabels
   colnames(responses) <- ilabels
    results <- list(scores=scores,responses = responses)} else {
    scores <- scoresList
     colnames(scores) <- ilabels
     results <- scores}
   #class(results) <- c("psych", "score.items")
    return(results)
 }

#created July 27, 2018
"scoreVeryFast" <- function(keys,items,totals=FALSE, min=NULL,max=NULL,count.responses=FALSE) {  #just scores by addition, no imputation, nothing fancy
if(is.data.frame(keys)) stop("I think you reversed keys and items.  I am stopping")
#use this for parallelism 
smallFunction <- function(scale,keys) {
          	pos.item <- items[,which(keys[,scale] > 0)]
            neg.item <- items[,which(keys[,scale] < 0)]
          	neg.item <- max + min - neg.item
           	sub.item <- cbind(pos.item,neg.item)
           if(count.responses) rs <- rowSums(!is.na(sub.item))
   if(totals) { scores <- rowSums(sub.item,na.rm=TRUE)} else {scores <- rowMeans(sub.item,na.rm=TRUE) } 
    if(count.responses) {return(c(scores,rs))} else {return(scores) }
  }
  
 if(is.list(keys)) {select <- sub("-","",unlist(keys))
   select <- select[!duplicated(select)]
   if(any( !(select %in% colnames(items)) )) {
 cat("\nVariable names are incorrect. Offending items are ", select[which(!(select %in% colnames(items)))],"\n")
 stop("Improper input.  See above. ")}
 } else {
      keys <- keys2list(keys)
     select <- selectFromKeyslist(colnames(items),keys)
      select <- select[!duplicated(select)]}
      items <- items[,select,drop=FALSE]
      n.subjects <- NROW(items)
      keys <- make.keys(items,keys)   #added 9/9/16 
   keys <- as.matrix(keys)   #just in case they were not matrices to start with
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
     abskeys <- abs(keys)
     keynames <- colnames(keys)
    if(is.null(keynames)) {if (totals) {keynames<- paste("S",1:n.keys,sep="")} else {
    	             keynames <- paste("A",1:n.keys,sep="")} }
    	             
   num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale


    n.subjects <- dim(items)[1]
    
   items <- as.matrix(items)
   scores <- matrix(NaN,ncol=n.keys,nrow=n.subjects)
  
    if (is.null(min)) {min <- min(items,na.rm=TRUE)}
    if (is.null(max)) {max <- max(items,na.rm=TRUE)}

     #use mapply for debugging, mcmapply for parallel processing
   #items are global and not passed
    scoresList <- mcmapply(smallFunction,c(1:n.keys),MoreArgs=list(keys=keys))  #the parallelized function
   
    if(count.responses) { scores <- scoresList[1:n.subjects,]
    responses <- scoresList[(n.subjects+1):nrow(scoresList),]
    
   colnames(scores) <- keynames
   colnames(responses) <- keynames
    results <- list(scores=scores,responses = responses)} else {
    
     scores <- scoresList
     colnames(scores) <- keynames
     results <- scores}
  
    return(results)
   }

 

