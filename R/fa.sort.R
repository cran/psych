#modified Sept 5, 2016 to sort Structure as well as loadings
"fa.sort" <- 
function(fa.results,polar=FALSE) {
  omega <- FALSE
  cor2 <- FALSE
  con.i <- FALSE
  fa.ci <- extension <-extend <- NA  #put in to avoid being identified as not defined.  Seems nuts
  Structure <- NULL  #in case we are not doing fa
#  if(length(class(fa.results)) > 1)  { value <- class(fa.results)[2] } else {value="other"}
 #This next section was added December 7, 2019 to change from class(x)[2] to inherits(x, ...)
 if(length(class(fa.results)) > 1)  {
    names <- cs(omega,omegaSem, fa.ci, iclust, fa, principal, extension, extend, cor2)
    value <- inherits(fa.results,names,which=TRUE)   # value <- class(x)[2]
    if(any(value > 1) ) { value <- names[which(value > 0)]} else {value <- "other"}
    
     } else {value <- "other"}
 

switch(value, 
omega =  { omega <- TRUE
          omegaSem <- FALSE
         factors <- as.matrix(fa.results$schmid$oblique)
          sl <- fa.results$schmid$sl},
        
omegaSem = {omega=TRUE
         omegaSem <- TRUE
         factors <- as.matrix(fa.results)
         sl <- as.matrix(fa.results)
         },

fa.ci = {factors <- fa.results$loadings
           if(!is.null(fa.results$Phi))  {Phi <-fa.results$Phi }
           con.i <-  TRUE
           ci <- fa.results$cis$ci
           cip <- fa.results$cis$p
           },
 
iclust = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
             
fa = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}
             Structure <- fa.results$Structure},
             
principal = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
             
extension = {factors <- as.matrix(fa.results$loadings)
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
             
extend = {factors <- as.matrix(fa.results$loadings)
               if(!is.null(fa.results$Structure)) Structure <- fa.results$Structure
             if(!is.null(fa.results$Phi))  {Phi <- fa.results$Phi}},
             
cor2 ={ factors <- as.matrix(fa.results$r)
        pval <- as.matrix(fa.results$pval) 
        cor2 <- TRUE
        },
             
other = {factors <- fa.results})

#now we have found the factor loadings  from the various possibilities
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
  		if(!is.null(Structure)) {Structure[1:nitems,] <- Structure[ord$ix,]
  		 rownames(Structure)[1:nitems] <- rownames(Structure)[ord$ix]
  		 }
  		
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
   				if(!is.null(Structure) ) {Structure[first:last,] <- Structure[item[ord$ix+first-1],]
   				rownames(Structure)[first:last] <- rownames(Structure)[ord$ix+first-1]
   				}
   				rownames(factors)[first:last] <- rownames(factors)[ord$ix+first-1]
   				
   			
   				if(con.i) { ci[first:last,] <- ci[item[ord$ix+first-1],] #if we are doing confidence intervals
                            cip[first:last,] <- cip[item[ord$ix+first-1],] }
        
   				total.ord[first:last] <- total.ord[ord$ix+first-1 ]
   		 		first <- first + items[i]  }
          		 
          		}  
 
 }
 
 if(cor2) {
          fa.results$r <- factors
           fa.results$pval <-  pval[total.ord,]
           }
         
 if(omega) {if(!omegaSem) fa.results$schmid$oblique <- factors
 
         #if the input was from omega, then sort the schmid leiman solution as well
         loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems))
         #nfactors <- dim(sl)[2]-4  #g, h2, u2, p2
       
        if("h2" %in% colnames(factors)) nfactors <- which(colnames(sl)=="h2") -1 
         if(omegaSem) nfactors <- NCOL(sl) -1
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
          		 
           if(omegaSem) {fa.results <- sl } else { fa.results$schmid$sl <- sl}  
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
         fa.results$uniquenesses <- fa.results$uniquenesses[total.ord]
         if(!is.null(Structure)) { fa.results$Structure <- Structure}
         }
          return(fa.results)
         }     
