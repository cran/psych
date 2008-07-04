"make.keys" <-
function(nvars,keys.list,key.labels=NULL,item.labels=NULL) {
 nkeys <- length(keys.list) 
keys <- matrix(rep(0,nvars*nkeys),ncol=nkeys)
for (i in 1:nkeys) {
		list.i <-  unlist(keys.list[[i]])
		keys[abs(list.i),i] <- sign(list.i ) 
			}
if(!is.null(key.labels)) {colnames(keys) <- key.labels} else {colnames(keys) <- names(keys.list)}
if(!is.null(item.labels)) {rownames(keys) <- item.labels} 
return(keys)}
#written June 11, 2008