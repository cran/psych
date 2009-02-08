"structure.list" <-
function(nvars,f.list,f=NULL,f.labels=NULL,item.labels=NULL) {
 nfactors <- length(f.list) 
fmodel <- matrix(rep(0,nvars*nfactors),ncol=nfactors)
for (i in 1:nfactors) {
        if(!is.null(f.list[[i]])) {
		list.i <-  unlist(f.list[[i]])
		fmodel[abs(list.i),i] <- paste(f,letters[i],list.i,sep="") } 
			}

if(!is.null(f.labels)) {colnames(fmodel) <- f.labels} else {colnames(fmodel) <- names(f.list)}
if(!is.null(item.labels)) {rownames(fmodel) <- item.labels} 
return(fmodel)}
#written Jan 22, 2009

"phi.list" <-
function(nf,f.list,f.labels=NULL) {
 nkeys <- length(f.list) 
phi <- diag(1,nf,nf)
for (i in 1:nkeys) {
		list.i <-  unlist(f.list[[i]])

		phi[list.i,i] <- paste("r",letters[i],letters[list.i],sep="") 
			}

if(!is.null(f.labels)) {colnames(phi) <- f.labels} else  {colnames(phi) <- paste("F",1:nf,sep="")}
rownames(phi) <- colnames(phi)
return(phi)}
#written Jan 22, 2009