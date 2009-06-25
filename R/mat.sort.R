"mat.sort" <-
function(m,f=NULL) {
if (is.null(f) ) {f <- fa(m) } 
if(is.list(f) && (!is.null(loadings(f)))) {load <- loadings(f)}  else {load <- f}
    load  <- as.matrix(load)
 	nitems <- NROW(load)
 	nfactors <-NCOL(load)
  	loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems),unclass(load))
#first sort them into clusters
      
  		#first find the maximum for each row and assign it to that cluster
  		 loads$cluster <- apply(abs(load),1,which.max)
 		 ord <- sort(loads$cluster,index.return=TRUE)
  		loads[1:nitems,] <- loads[ord$ix,]
 		rownames(loads)[1:nitems] <- rownames(loads)[ord$ix]
 		 
  #now sort column wise
  #now sort the loadings that have their highest loading on each cluster
  
  		items <- table(loads$cluster)   #how many items are in each cluster?
  		first <- 1
  		item <- loads$item
		for (i in 1:length(items)) {
		if(items[i] > 0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(loads[first:last,i+2]),decreasing=TRUE,index.return=TRUE)
   				loads[first:last,3:(nfactors+2)] <- load[item[ord$ix+first-1],]
   				loads[first:last,1] <- item[ord$ix+first-1]
   				rownames(loads)[first:last] <- rownames(loads)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
      item.order <- loads[,1]    		 
   m  <- m[item.order,item.order]
   return(m)
   }
	