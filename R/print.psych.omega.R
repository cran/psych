"print.psych.omega" <- 
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

if(is.null(cut)) cut <- .2
	 cat( x$title,"\n") 
	 cat("Call: ")
     print(x$call)
 	cat("Alpha:                ",round(x$alpha,digits),"\n") 
 	cat("G.6:                  ",round(x$G6,digits),"\n")
 	cat("Omega Hierarchical:   " ,round(x$omega_h,digits),"\n")
 	cat("Omega H asymptotic:   " ,round(x$omega.lim,digits),"\n")
 	cat("Omega Total           " ,round(x$omega.tot,digits),"\n")
            
	cat("\nSchmid Leiman Factor loadings greater than ",cut, "\n")
	
	   loads <- x$schmid$sl
	   nfactor <- ncol(loads)-2
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
	        h2 <- round(rowSums(load.2^2),digits)
	        loads <-round(loads,digits)
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")  
         	fx[abs(loads)< cut] <- paste(rep(" ", nc), collapse = "")
            h2 <- round(rowSums(load.2^2),digits)
         	u2 <- 1 - h2
	    	print(cbind(fx[,1:(nfactor+sort)],h2,u2),quote="FALSE")
	    	#print(fx,quote="FALSE")
	   
       numfactors <- dim(x$schmid$sl)[2] -2
       eigenvalues <- diag(t(x$schmid$sl[,1:numfactors]) %*% x$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       ev.rnd <- round(eigenvalues,digits)
       print(ev.rnd,digits=digits)
      
  	 maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
  	 gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
  	
   	cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits))
   	cat("\nThe degrees of freedom for the model is",x$schmid$dof," and the fit was ",round(x$schmid$objective,digits),"\n")
   	if(!is.na(x$schmid$n.obs)) {cat("The number of observations was ",x$schmid$n.obs, " with Chi Square = ",round(x$schmid$STATISTIC,digits), " with prob < ", signif(x$schmid$PVAL,digits),"\n")}
  stats.df <- t(data.frame(sqrt(x$stats$R2),x$stats$R2,2*x$stats$R2 -1))
  cat("\nMeasures of factor score adequacy             \n")
  rownames(stats.df) <- c("Correlation of scores with factors  ","Multiple R square of scores with factors ","Minimum correlation of factor score estimates")
  print(round(stats.df,digits))
  #cat("\nMeasures of factor score adequacy             ",names(eigenvalues))
  #cat("\nCorrelation of scores with factors            ",round(sqrt(x$stats$R2),digits))
  #cat("\nMultiple R square of scores with factors      " ,round(x$stats$R2,digits))
  #cat("\nMinimum correlation of factor score estimates ", round(2*x$stats$R2 -1,digits),"\n")
   
   }