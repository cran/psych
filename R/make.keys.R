"make.keys" <-
function(nvars,keys.list,item.labels=NULL,key.labels=NULL) {
if(!is.null(ncol(nvars)))  {item.labels <- colnames(nvars)
                         nvars <-  ncol(nvars)} else {
if(!is.numeric(nvars)) {item.labels <- nvars 
                      nvars <- length(item.labels)} }
nkeys <- length(keys.list) 
keys <- matrix(rep(0,nvars*nkeys),ncol=nkeys)
for (i in 1:nkeys) {
		list.i <-  unlist(keys.list[[i]])
		if((is.character(list.i)) && !is.null(item.labels)) {
		neg <- grep("-",list.i)
		list.i <- sub("-","",list.i)
		list.i <- match(list.i,item.labels)
		if(!any(is.na(neg))) list.i[neg] <- -list.i[neg]}
		keys[abs(list.i),i] <- sign(list.i ) 
			}
if(!is.null(key.labels)) {colnames(keys) <- key.labels} else {colnames(keys) <- names(keys.list)}
if(!is.null(item.labels)) {rownames(keys) <- item.labels} 
return(keys)}
#written June 11, 2008
#revised Sept 15, 2013 to allow for symbolic keys





