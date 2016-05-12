"iclust.sort"<- function (ic.load,cut=0,labels=NULL,keys=FALSE, clustsort=TRUE) {ICLUST.sort(ic.load,labels,keys,clustsort)}
"ICLUST.sort"<- function (ic.load,cut=0,labels=NULL,keys=FALSE, clustsort=TRUE) {
     if(is.matrix(ic.load)) {loadings <- ic.load
     pattern <- as.matrix(loadings)} else { loadings <- ic.load$loadings
                                 pattern <- as.matrix(ic.load$pattern)}
     
     nclust <- dim(loadings)[2]
     nitems <- dim(loadings)[1]
     loadings <- as.matrix(loadings)   #just in case there is just one cluster
     loadings <- unclass(loadings)     #to get around the problem of a real loading matrix
     if(nclust > 1) {eigenvalue <- diag(t(pattern)%*% loadings)  #put the clusters into descending order by eigenvalue
                     evorder <- order(eigenvalue,decreasing=TRUE)
                     if(clustsort) loadings <- loadings[,evorder]  #added the clustsort option 2011.12.22 until now had always sorted
                     }
    if (length(labels)==0) {
    var.labels <- rownames(loadings)} else {var.labels=labels}
    if (length(var.labels)==0) {var.labels =paste('V',seq(1:nitems),sep='')} #unlabled variables

   
  loads <- data.frame(item=seq(1:nitems),content=var.labels,cluster=rep(0,nitems),loadings)
 
 
  #first find the maximum for each row and assign it to that cluster
   loads$cluster <- apply(abs(loadings),1,which.max)
  for (i in 1:nitems) {if (abs(loadings[i,loads$cluster[i]]) < cut) {loads$cluster[i] <- nclust+1}} #assign the ones that missed the cut a location
 
  ord <- sort(loads$cluster,index.return=TRUE)
  loads[1:nitems,] <- loads[ord$ix,]
  rownames(loads)[1:nitems] <- rownames(loads)[ord$ix]
  
  items <- c(table(loads$cluster),1)   #how many items are in each cluster?
  if(length(items) < (nclust+1)) {items <- rep(0,(nclust+1))   #this is a rare case where some clusters don't have anything in them
    for (i in 1:nclust+1) {items[i] <- sum(loads$cluster==i) }  }

  #now sort the loadings that have their highest loading on each cluster
   first <- 1
	for (i in 1:nclust) {
	if(items[i]>0 ) {
	last <- first + items[i]- 1
	ord <- sort(abs(loads[first:last,i+3]),decreasing=TRUE,index.return=TRUE)
   loads[first:last,] <- loads[ord$ix+first-1,]
    rownames(loads)[first:last] <- rownames(loads)[ord$ix+first-1]
    first <- first + items[i]}
    }
    if (first < nitems) loads[first:nitems,"cluster"] <- 0   #assign items less than cut to 0
      if(keys) {result <- list(sorted=loads,clusters=factor2cluster(loadings))} else  result <- list(sorted=loads)
   class(result) <- c("psych","iclust.sort")   #need to clean up print to make this work
   return(result)
}
 #revised August 8, 2007 to add cluster keying option and to allow us to work with factor analysis output
 #revised Sept 15, 2007 to remove the "loadings" parameter
 #revised Ausgust 30, 2008 to make class psych
 #revised August 28, 2012 to meet a request from Gundmundur Arnkelsson to be able to print from principal output.
 