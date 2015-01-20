"print.psych.schmid" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

if(is.null(cut)) cut <- .2
	 cat( "Schmid-Leiman analysis \n") 
	 cat("Call: ")
     print(x$Call)
 
	cat("\nSchmid Leiman Factor loadings greater than ",cut, "\n")
	
	   loads <- x$sl
	   nfactor <- ncol(loads)-3
	   if(sort) {
	      ord <- sort(abs(loads[,1]),decreasing=TRUE,index.return=TRUE)
	      loads[,] <- loads[ord$ix,]
	      rownames(loads) <- rownames(loads)[ord$ix]
	      loads <- cbind(v=ord$ix,loads)
	   } #end sort 
	        tn <- colnames(loads)
	        loads <- data.frame(loads)
	        colnames(loads) <- tn  #this seems weird, but otherwise we lose the F* name
	       
	        if(sort) {loads[,1] <- as.integer(loads[,1])
	                 load.2 <- loads[,2:(nfactor+1)]} else {load.2 <- loads[,1:nfactor] }    
	        h2 <- round(loads[,"h2"],digits)
	        u2 <- round(loads[,"u2"],digits) 
	        loads <- round(loads,digits)
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")  
         	fx[abs(loads)< cut] <- paste(rep(" ", nc), collapse = "")
           # h2 <- round(rowSums(load.2^2),digits)
           #u2 <- 1 - h2
         	
         	p2 <- loads[,"p2"]
         	mp2 <- mean(p2)
         	vp2 <- var(p2)
         	p2 <- round(p2,digits)
	    	print(cbind(fx[,1:(nfactor+sort)],h2,u2,p2),quote="FALSE")
	    	
	   
       numfactors <- dim(x$sl)[2] -3
       eigenvalues <- diag(t(x$sl[,1:numfactors]) %*% x$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       ev.rnd <- round(eigenvalues,digits)
       print(ev.rnd,digits=digits)
      
  	 maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
  	 gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
  	 
  	 
  		cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits))
    cat("\nmean percent general = ",round(mp2,digits), "   with sd = ", round(sqrt(vp2),digits), "and cv of ",round(sqrt(vp2)/mp2,digits),"\n")
   
cat("\n The orthogonal loadings were \n")

	load <- x$orthog   
cut <- 0
 
    #but, if we are print factors of covariance matrices, they might be very small
     #       cut <- min(cut,max(abs(load))/2)   #removed following a request by  Reinhold Hatzinger
 	nitems <- dim(load)[1]
 	nfactors <- dim(load)[2]
  	loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems),unclass(load))
  	u2.order <- 1:nitems  #used if items are sorted
 if(sort) {
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
  		
		for (i in 1:length(items)) {# i is the factor number
		if(items[i] > 0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(loads[first:last,i+2]),decreasing=TRUE,index.return=TRUE)
				u2.order[first:last] <- item[ord$ix+first-1]
   				loads[first:last,3:(nfactors+2)] <- load[item[ord$ix+first-1],]
   				loads[first:last,1] <- item[ord$ix+first-1]
   				rownames(loads)[first:last] <- rownames(loads)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
         }    #end of sort 		 
    #they are now sorted, don't print the small loadings if cut > 0 

          	ncol <- dim(loads)[2]-2
          	rloads <- round(loads,digits)
	    	fx <- format(rloads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")
	    	 fx.1 <- fx[,1,drop=FALSE]    #drop = FALSE  preserves the rownames for single factors
	    	 fx.2 <- fx[,3:(2+ncol)]
	    	 load.2 <- as.matrix(loads[,3:(ncol+2)])
         	fx.2[abs(load.2) < cut] <- paste(rep(" ", nc), collapse = "")
         if(sort) {	fx <- data.frame(V=fx.1,fx.2)
         	if(dim(fx)[2] <3) colnames(fx) <- c("V",colnames(x$loadings)) #for the case of one factor
         } else {fx <- data.frame(fx.2)
            colnames(fx) <- colnames(x$orthog)}
         	if(nfactors > 1) {if(is.null(x$Phi)) {h2 <- rowSums(load.2^2)} else {h2 <- diag(load.2 %*% x$Phi %*% t(load.2)) }} else {h2 <-load.2^2}
         	if(!is.null(x$uniquenesses)) {u2 <- x$uniquenesses[u2.order]}  else {u2 <- u2[u2.order]}
         	#h2 <- round(h2,digits)
         	vtotal <- sum(h2 + u2)
           if(isTRUE(all.equal(vtotal,nitems))) {
           cat("Standardized loadings based upon correlation matrix\n")
	    	print(cbind(fx,h2,u2),quote="FALSE",digits=digits) } else {
	    	cat("Unstandardized loadings based upon covariance matrix\n") 
	    	print(cbind(fx,h2,u2,H2=h2/(h2+u2),U2=u2/(h2+u2)),quote="FALSE",digits=digits)}
	   
 		
      	   #adapted from print.loadings
      	  if(is.null(x$Phi)) {if(nfactors > 1)  {vx <- colSums(load.2^2) } else {vx <- sum(load.2^2)
      	                                             }} else {vx <- diag(x$Phi %*% t(load) %*% load)
      	                                                      }
          
          names(vx) <- colnames(x$orthog)
          varex <- rbind("SS loadings" =   vx)
          varex <- rbind(varex, "Proportion Var" =  vx/vtotal)
         if (nfactors > 1) 
                      varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/vtotal))
                      
             cat("\n")
  
    print(round(varex, digits))


   	
   	 if(!is.null(x$dof)) {cat("\nThe degrees of freedom are",x$dof," and the fit is ",round(x$objective,digits),"\n")
    if(!is.null(x$n.obs)&&!is.na(x$n.obs)) {cat("The number of observations was ",x$n.obs, " with Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits))}
    }
    if(!is.null(x$rms)) {cat("\nThe root mean square of the residuals is ", round(x$rms,digits),"\n") }
     if(!is.null(x$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$crms,digits)) }
     if(!is.null(x$RMSEA)) {cat("\nRMSEA index = ",round(x$RMSEA[1],digits+1), " and the", (1- x$RMSEA[4])*100,"% confidence intervals are ",round(x$RMSEA[2:3],digits+1))  }
   	if(!is.null(x$BIC)) {cat("\nBIC = ",round(x$BIC,digits))}
   	

   }
