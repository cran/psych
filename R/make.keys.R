"make.keys" <-
function(nvars,keys.list,item.labels=NULL,key.labels=NULL) {
if(!is.null(ncol(nvars)))  {item.labels <- colnames(nvars)
                         nvars <-  ncol(nvars)} else {
if(!is.numeric(nvars)) {item.labels <- nvars 
                      nvars <- length(item.labels)} }
nkeys <- length(keys.list) 
keys <- matrix(rep(0,nvars*nkeys),ncol=nkeys)
for (i in 1:nkeys) {
      if(!is.null(keys.list[[i]])) {
		list.i <-  unlist(keys.list[[i]])
		if((is.character(list.i)) && !is.null(item.labels)) {
		neg <- grep("-",list.i)
		list.i <- sub("-","",list.i)
		list.i <- match(list.i,item.labels)
		if(!any(is.na(neg))) list.i[neg] <- -list.i[neg]}
		keys[abs(list.i),i] <- sign(list.i ) 
			}
			}
if(!is.null(key.labels)) {colnames(keys) <- key.labels} else {colnames(keys) <- names(keys.list)}
if(!is.null(item.labels)) {rownames(keys) <- item.labels} 
return(keys)}
#written June 11, 2008
#revised Sept 15, 2013 to allow for symbolic keys
#revised November 21, 2018 to allow null keys



# 
#    "keys2list" <- function(keys,sign=TRUE) {
#       keys.list <- list()
#       nkeys <- ncol(keys)
#       for (i in 1:nkeys) {temp <- rownames(keys)[which(keys[,i] < 0)]
#      if(sign && (length(temp)  >0)) temp <- paste0("-",temp)
#       keys.list[[i]] <- c(rownames(keys)[which(keys[,i] > 0)],temp) 
#       } 
#       names(keys.list) <- colnames(keys)
#       keys.list}
      
#Added July 9, 2017
"selectFromKeys" <- function(keys.list) {
    select <- sub("-","",unlist(keys.list))
      select <- select[!duplicated(select)]
      return(select) 
      }

#Basically, the opposite of make.keys
#Takes a keys matrix and converts it to a list structure (with negative signs appropriately placed)   
#9/10/16
#revised 6/10/18 to not change the order of keys

   "keys2list" <- function(keys,sign=TRUE) {
      keys.list <- list()
      nkeys <- ncol(keys)
      for (i in 1:nkeys) {temp <- rownames(keys)[which(keys[,i] != 0)]
      wk <- which(keys[,i] < 0)
      temp[temp %in% names(wk)] <- paste0("-",temp[temp %in% names(wk)])
      keys.list[[i]] <- temp
      
     #if(sign && (length(temp)  >0)) temp <- paste0("-",temp)
     # keys.list[[i]] <- c(rownames(keys)[which(keys[,i] > 0)],temp) 
      } 
      names(keys.list) <- colnames(keys)
      keys.list}
      
    
"makePositiveKeys" <- function(keys.list,sign=FALSE) {
     nkeys <- length(keys.list)
     pos.list <- neg.list  <- list()
     for(i in 1:nkeys) {temp <- keys.list[[i]]
      neg.pos <- grep("-",temp)
      if(length(neg.pos)  >0) {
      pos.list[[i]]  <- temp[-neg.pos]
      neg.list[[i]]  <- temp[neg.pos]} else {pos.list[[i]] <- keys.list[[i]]
      neg.list[[i]] <-keys.list[[i]]}
      if(!sign) neg.list[[i]] <- sub("-","",neg.list[[i]])
      }
    names(pos.list) <- names(keys.list)
    names(neg.list) <- paste0(names(keys.list),".neg")
    return(c(pos.list,neg.list)) 
    }
      
    
   