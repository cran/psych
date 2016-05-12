"fa.sort" <- 
function(fa.results,polar=FALSE) {
  omega <- FALSE
  con.i <- FALSE
  if(length(class(fa.results)) > 1)  { value <- class(fa.results)[2] } else {value="other"}

switch(value, 
omega =  { omega <- TRUE
        factors <- as.matrix(fa.results$schmid$oblique)
        sl <- fa.results$schmid$sl},

fa.ci = {factors <- fa.results$loadings
           if(!is.null(fa.results$Phi))  {Phi <-fa.results$Phi }
           con.i <-  TRUE
           ci <- fa.results$cis$ci
           cip <- fa.results$cis$p
           },
 
iclust = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
fa = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
principal = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
extension = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
extend = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
other = {factors <- fa.results})

nitems <- dim(factors)[1]
nfactors <- dim(factors)[2]
total.ord <- rep(NA,nitems)
if(polar) { pol.ord <- polar(factors)[,1]
      factors[1:nitems,] <- factors[pol.ord,]
 		rownames(factors)[1:nitems] <- rownames(factors)[pol.ord ] } else {
 		
if(is.null(rownames(factors))) {rownames(factors) <- paste("V",1:nitems)}

loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems))

 		#first sort them into clusters
  		#first find the maximum for each row and assign it to that cluster
  		 loads$cluster <- apply(abs(factors),1,which.max)
 		 ord <- sort(loads$cluster,index.return=TRUE)
  		factors[1:nitems,] <- factors[ord$ix,]
 		rownames(factors)[1:nitems] <- rownames(factors)[ord$ix]
        total.ord <- ord$ix
        if(con.i) { ci[1:nitems,] <- ci[ord$ix,] #if we are doing confidence intervals
                                   cip[1:nitems,] <- cip[ord$ix,] }
 		
 		 
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
   			
   				if(con.i) { ci[first:last,] <- ci[item[ord$ix+first-1],] #if we are doing confidence intervals
                            cip[first:last,] <- cip[item[ord$ix+first-1],] }
        
   				total.ord[first:last] <- total.ord[ord$ix+first-1 ]
   		 		first <- first + items[i]  }
          		 
          		}  
 
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
              
         } else {
            if((!is.matrix(fa.results)) && (!is.data.frame(fa.results)))  {fa.results$loadings <- factors
             if(con.i) { rownames(ci) <- rownames(factors)
                        fa.results$ci <- ci
                      rownames(cip) <- rownames(factors)
                      colnames(cip) <- colnames(factors)
                      fa.results$cip <- cip}
                       
              } else {
              fa.results <- factors} }
        #note that h2 and complexities were not sorted, we need to do this now
       if(is.list(fa.results))  {fa.results$order <- total.ord 
         fa.results$complexity <- fa.results$complexity[total.ord]
         fa.results$communality <- fa.results$communality[total.ord]
         fa.results$uniquenesses <- fa.results$uniquenesses[total.ord]}
          return(fa.results)
         }     
