"print.factor.pa" <-
function(x,digits=2,all=FALSE,cutoff=NULL,sort=FALSE,...) { 

 if(is.null(cutoff)) cutoff <- .3
 	load <- x$loadings
 	nitems <- dim(load)[1]
 	nfactors <- dim(load)[2]
  	loads <- data.frame(item=seq(1:nitems),cluster=rep(0,nitems),unclass(load))
 if(sort) {
 		#first sort them into clusters
  		#first find the maximum for each row and assign it to that cluster
  		 loads$cluster <- apply(abs(load),1,which.max)
 		 ord <- sort(loads$cluster,index.return=TRUE)
  		loads[1:nitems,] <- loads[ord$ix,]
 		rownames(loads)[1:nitems] <- rownames(loads)[ord$ix]
 		 
  #now sort column wise
   		items <- c(table(loads$cluster),1)   #how many items are in each cluster?
  		if(length(items) < (nfactors+1)) {items <- rep(0,(nfactors+1))   #this is a rare case where some clusters don't have anything in them
   		                                  for (i in 1:nfactors+1) {items[i] <- sum(loads$cluster==i) }  }

  #now sort the loadings that have their highest loading on each cluster
  		first <- 1
		for (i in 1:nfactors) {
		if(items[i]>0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(loads[first:last,i+2]),decreasing=TRUE,index.return=TRUE)
   				loads[first:last,] <- loads[ord$ix+first-1,]
   				 rownames(loads)[first:last] <- rownames(loads)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
         }    #end of sort 		 
    #they are now sorted, don't print the small loadings 
          	ncol <- dim(loads)[2]-2
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")
	    	 fx.1 <- fx[,1]
	    	 fx.2 <- fx[,3:(2+ncol)]
	    	 load.2 <- loads[,3:(ncol+2)]
         	fx.2[abs(load.2)< cutoff] <- paste(rep(" ", nc), collapse = "")
         	fx <- data.frame(V=fx.1,fx.2)
	    	print(fx,quote="FALSE")
 		
      	   #adapted from print.loadings
      	   vx <- colSums(load.2^2)
           varex <- rbind("SS loadings" =   vx)
            varex <- rbind(varex, "Proportion Var" =  vx/nitems)
            if (nfactors > 1) 
                      varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/nitems))
  
    cat("\n")
  
    print(round(varex, digits))
    
    if(!is.null(x$phi))  { 
       cat ("\n With factor correlations of \n" )
       colnames(x$phi) <- rownames(x$phi) <- colnames(x$loadings)
       print(round(x$phi,digits))} else {
       if(!is.null(x$rotmat)) {
             U <- x$rotmat
           phi <- t(U) %*% U
          phi <- cov2cor(phi) 
           cat ("\n With factor correlations of \n" )
       colnames(phi) <- rownames(phi) <- colnames(x$loadings)
       print(round(phi,digits))
            } }
            
       objective <- x$criteria[1]
     if(!is.null(objective)) {    cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.\n")
    cat("\nThe degrees of freedom for the model is",x$dof," and the fit was ",round(objective,digits),"\n") 
   	if(!is.na(x$n.obs)) {cat("The number of observations was ",x$n.obs, " with Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits),"\n")}

}
          
	 }  
