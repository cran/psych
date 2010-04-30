"fa.sort" <- 
function(fa.results) {
  omega <- FALSE
if(!is.null(class(fa.results))) {if(class(fa.results)[2] =="omega") {
 omega <- TRUE
factors <- as.matrix(fa.results$schmid$oblique)
sl <- fa.results$schmid$sl
} else {  #the normal factor case

 Phi <- NULL  #the default case
 if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {factors <- as.matrix(fa.results$loadings)
                if(!is.null(fa.results$Phi)) Phi <- fa.results$Phi} else {factors <- fa.results}
 }
}                
nitems <- dim(factors)[1]
nfactors <- dim(factors)[2]
if(is.null(rownames(factors))) {rownames(factors) <- paste("V",1:nitems)}

loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems))

 		#first sort them into clusters
  		#first find the maximum for each row and assign it to that cluster
  		 loads$cluster <- apply(abs(factors),1,which.max)
 		 ord <- sort(loads$cluster,index.return=TRUE)
  		factors[1:nitems,] <- factors[ord$ix,]
 		rownames(factors)[1:nitems] <- rownames(factors)[ord$ix]
 		
 		 
  #now sort column wise
  #now sort the loadings that have their highest loading on each cluster
  
  		items <- table(loads$cluster)   #how many items are in each cluster?
  		first <- 1
  		item <- loads$item
		for (i in 1:length(items)) {
		if(items[i] > 0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(factors[first:last,i]),decreasing=TRUE,index.return=TRUE)
   				factors[first:last,] <- factors[item[ord$ix+first-1],]
   				loads[first:last,1] <- item[ord$ix+first-1]
   				rownames(factors)[first:last] <- rownames(factors)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
          		 
         
 if(omega) {fa.results$schmid$oblique <- factors
 
         #if the input was from omega, then sort the schmid leiman solution as well
         loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems))
         nfactors <- dim(sl)[2]-4  #g, h2, u2, p2
         if(nfactors > 1) {
         loads$cluster <- apply(abs(sl[,2:(nfactors+1)]),1,which.max) +1} else {loads$cluster <- rep(1,nitems) }
 		 ord <- sort(loads$cluster,index.return=TRUE)
  		sl[1:nitems,] <- sl[ord$ix,]
 		rownames(sl)[1:nitems] <- rownames(sl)[ord$ix]
 		items <- table(loads$cluster)   #how many items are in each cluster?
  		first <- 1
  		item <- loads$item
		for (i in 1:length(items)) {
		if(items[i] > 0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(sl[first:last,i+1]),decreasing=TRUE,index.return=TRUE)
   				sl[first:last,] <- sl[item[ord$ix+first-1],]
   				loads[first:last,1] <- item[ord$ix+first-1]
   				rownames(sl)[first:last] <- rownames(sl)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
          		 
          fa.results$schmid$sl <- sl
         
         
         
         
         
         
         
         
         
         
         } else {if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {fa.results$loadings <- factors} else {
              fa.results <- factors} }
              
          return(fa.results)
         }     
