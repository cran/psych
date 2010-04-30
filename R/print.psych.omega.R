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
	        h2 <- round(rowSums(load.2^2),digits)
	        loads <-round(loads,digits)
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")  
         	fx[abs(loads)< cut] <- paste(rep(" ", nc), collapse = "")
            h2 <- round(rowSums(load.2^2),digits)
         	u2 <- 1 - h2
         	
         	p2 <- loads[,ncol(loads)]
         	mp2 <- mean(p2)
         	vp2 <- var(p2)
         	p2 <- round(p2,digits)
	    	print(cbind(fx[,1:(nfactor+sort)],h2,u2,p2),quote="FALSE")
	    	#print(fx,quote="FALSE")
	   
       numfactors <- dim(x$schmid$sl)[2] -3
       eigenvalues <- diag(t(x$schmid$sl[,1:numfactors]) %*% x$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       ev.rnd <- round(eigenvalues,digits)
       print(ev.rnd,digits=digits)
      
  	 maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
  	 gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
  	
   	cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits))
    cat("\nmean percent general = ",round(mp2,digits), "   with sd = ", round(sqrt(vp2),digits), "and cv of ",round(sqrt(vp2)/mp2,digits),"\n")
   
   if(FALSE) {cat("\nCompare three different sets of goodness of fit measures.\nFirst, the omega model:")   #this next section needs to be rethought 
   
   if(!is.null(x$stats$dof)) {cat("\nThe degrees of freedom for the omega model are",x$stats$dof," and the fit was ",round(x$stats$objective,digits),"\n")
   if(!is.null(x$stats$n.obs)&&!is.na(x$stats$n.obs)) {cat("The number of observations was ",x$stats$n.obs, " with Chi Square = ",round(x$stats$STATISTIC,digits), " with prob < ", signif(x$schmid$PVAL,digits),"\n")}
    }
    if(!is.null(x$stats$rms)) {cat("\nThe root mean square of the residuals is ", round(x$stats$rms,digits),"\n") }
     if(!is.null(x$stats$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$stats$crms,digits),"\n") }
    if(!is.null(x$stats$RMSEA)) {cat("\nRMSEA and the ",x$stats$RMSEA[4]  ,"confidence intervals are ",round(x$stats$RMSEA[1:3],digits+1))  }
   	if(!is.null(x$stats$BIC)) {cat("\nBIC = ",round(x$stats$BIC,digits))}
   	cat("\nThen the simple factor model with no higher order factor")
   	}   #end commented out section
   	
   	 if(!is.null(x$schmid$dof)) {cat("\nThe degrees of freedom are",x$schmid$dof," and the fit is ",round(x$schmid$objective,digits),"\n")
    if(!is.null(x$schmid$n.obs)&&!is.na(x$schmid$n.obs)) {cat("The number of observations was ",x$schmid$n.obs, " with Chi Square = ",round(x$schmid$STATISTIC,digits), " with prob < ", signif(x$schmid$PVAL,digits))}
    }
    if(!is.null(x$schmid$rms)) {cat("\nThe root mean square of the residuals is ", round(x$schmid$rms,digits),"\n") }
     if(!is.null(x$schmid$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$schmid$crms,digits)) }
    if(!is.null(x$schmid$RMSEA)) {cat("\nRMSEA and the ",x$stats$RMSEA[4]  ,"confidence intervals are ",round(x$schmid$RMSEA[1:3],digits+1))  }
   	if(!is.null(x$schmid$BIC)) {cat("\nBIC = ",round(x$schmid$BIC,digits))}
   	
   	cat("\n\nCompare this with the adquacey of just a general factor and no group factors")
   	if(!is.null(x$gstats$dof)) {cat("\nThe degrees of freedom for just the general factor are",x$gstats$dof," and the fit is ",round(x$gstats$objective,digits),"\n")
    if(!is.null(x$gstats$n.obs)&&!is.na(x$gstats$n.obs)) {cat("The number of observations was ",x$gstats$n.obs, " with Chi Square = ",round(x$gstats$STATISTIC,digits), " with prob < ", signif(x$gstats$PVAL,digits))}
    }
    if(!is.null(x$gstats$rms)) {cat("\nThe root mean square of the residuals is ", round(x$gstats$rms,digits),"\n") }
     if(!is.null(x$gstats$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$gstats$crms,digits),"\n") }
    if(!is.null(x$gstats$RMSEA)) {cat("\nRMSEA and the ",x$gstats$RMSEA[4]  ,"confidence intervals are ",round(x$gstats$RMSEA[1:3],digits+1))  }
   	if(!is.null(x$gstats$BIC)) {cat("\nBIC = ",round(x$gstats$BIC,digits),"\n")}
   	
   
   	
   stats.df <- t(data.frame(sqrt(x$stats$R2),x$stats$R2,2*x$stats$R2 -1)) 
   
   
          
  cat("\nMeasures of factor score adequacy             \n")
  rownames(stats.df) <- c("Correlation of scores with factors  ","Multiple R square of scores with factors ","Minimum correlation of factor score estimates")
  print(round(stats.df,digits))
  #cat("\nMeasures of factor score adequacy             ",names(eigenvalues))
  #cat("\nCorrelation of scores with factors            ",round(sqrt(x$stats$R2),digits))
  #cat("\nMultiple R square of scores with factors      " ,round(x$stats$R2,digits))
  #cat("\nMinimum correlation of factor score estimates ", round(2*x$stats$R2 -1,digits),"\n")
   
   }