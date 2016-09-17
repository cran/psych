#created 7/12/16 
#just score items without a lot of stats
#basically scoreItems with all  the stats removed
"scoreFast"  <-
 function (keys,items,totals=FALSE,ilabels=NULL, missing=TRUE, impute="none",delete=TRUE,  min=NULL,max=NULL,digits=2) {
   cl <- match.call()
   raw.data <- TRUE
   if(impute == FALSE)  impute <- "none"
   if(is.list(keys)) {select <- sub("-","",unlist(keys))
      select <- select[!duplicated(select)]
      items <- items[select]
      keys <- make.keys(items,keys)}   #added 9/9/16 
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
        	#C <- cov(items,use="pairwise")
         # cov.scales  <- cov(scores,use="pairwise")    #and total scale variance
         # cov.scales2 <- diag(t(abskeys) %*% C^2 %*% abskeys)   # sum(C^2)  for finding ase
        }  else { #handle the case of missing data without imputation
           scores <- matrix(NaN,ncol=n.keys,nrow=n.subjects)
           totals <- FALSE  #just in case it was not already false
           #we could try to parallelize this next loop
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
           # C <- cov(items,use="pairwise")
            #cov.scales <- t(keys) %*% C %*% keys
           # cov.scales2 <- diag(t(abskeys) %*% C^2 %*% abskeys)  # sum(C^2)  for finding ase
          #  raw.data <- FALSE
         }  #end of treating missing without imputation


                   
    slabels <- colnames(keys)
    if (is.null(slabels)) {
    	if (totals) {slabels<- paste("S",1:n.keys,sep="")} else {
    	             slabels <- paste("A",1:n.keys,sep="")} }
   
   colnames(scores) <- slabels
    
   
  results <- scores
   #class(results) <- c("psych", "score.items")
    return(results)
 }
