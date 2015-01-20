"print.psych.omega" <- 
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 
xx <- x
if(!is.null(x$ci)) {
x <- x$om} 
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
	         h2 <- round(loads[,"h2"],digits)
	         u2 <- round(loads[,"u2"],digits) 
	        loads <-round(loads,digits)
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")  
         	fx[abs(loads)< cut] <- paste(rep(" ", nc), collapse = "")
           
         	
            p2 <- loads[,"p2"]
         	mp2 <- mean(p2)
         	vp2 <- var(p2)
         	p2 <- round(p2,digits)
	    	print(cbind(fx[,1:(nfactor+sort)],h2,u2,p2),quote="FALSE")
	    	
	   
       numfactors <- dim(x$schmid$sl)[2] -3
       eigenvalues <- diag(t(x$schmid$sl[,1:numfactors]) %*% x$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       ev.rnd <- round(eigenvalues,digits)
       print(ev.rnd,digits=digits)
      
  	 maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
  	 gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
  	
   	cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits))
    cat("\nmean percent general = ",round(mp2,digits), "   with sd = ", round(sqrt(vp2),digits), "and cv of ",round(sqrt(vp2)/mp2,digits),"\n")
   if(!is.null(x$ECV))  cat("Explained Common Variance of the general factor = ", round(x$ECV,digits),"\n")
   

   	
   	 if(!is.null(x$schmid$dof)) {cat("\nThe degrees of freedom are",x$schmid$dof," and the fit is ",round(x$schmid$objective,digits),"\n")
    if(!is.null(x$schmid$n.obs)&&!is.na(x$schmid$n.obs)) {cat("The number of observations was ",x$schmid$n.obs, " with Chi Square = ",round(x$schmid$STATISTIC,digits), " with prob < ", signif(x$schmid$PVAL,digits))}
    }
    if(!is.null(x$schmid$rms)) {cat("\nThe root mean square of the residuals is ", round(x$schmid$rms,digits),"\n") }
     if(!is.null(x$schmid$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$schmid$crms,digits)) }
     if(!is.null(x$schmid$RMSEA)) {cat("\nRMSEA index = ",round(x$schmid$RMSEA[1],digits+1), " and the", (1- x$schmid$RMSEA[4])*100,"% confidence intervals are ",round(x$schmid$RMSEA[2:3],digits+1))  }
   	if(!is.null(x$schmid$BIC)) {cat("\nBIC = ",round(x$schmid$BIC,digits))}
   	
   	cat("\n\nCompare this with the adequacy of just a general factor and no group factors")
   	if(!is.null(x$gstats$dof)) {cat("\nThe degrees of freedom for just the general factor are",x$gstats$dof," and the fit is ",round(x$gstats$objective,digits),"\n")
    if(!is.null(x$gstats$n.obs)&&!is.na(x$gstats$n.obs)) {cat("The number of observations was ",x$gstats$n.obs, " with Chi Square = ",round(x$gstats$STATISTIC,digits), " with prob < ", signif(x$gstats$PVAL,digits))}
    }
    if(!is.null(x$gstats$rms)) {cat("\nThe root mean square of the residuals is ", round(x$gstats$rms,digits),"\n") }
     if(!is.null(x$gstats$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$gstats$crms,digits),"\n") }
     if(!is.null(x$gstats$RMSEA)) {cat("\nRMSEA index = ",round(x$gstats$RMSEA[1],digits+1), " and the", (1- x$gstats$RMSEA[4])*100,"% confidence intervals are ",round(x$gstats$RMSEA[2:3],digits+1))  }
   	if(!is.null(x$gstats$BIC)) {cat("\nBIC = ",round(x$gstats$BIC,digits),"\n")}
   	
   
   	
   stats.df <- t(data.frame(sqrt(x$stats$R2),x$stats$R2,2*x$stats$R2 -1)) 
   
   
          
  cat("\nMeasures of factor score adequacy             \n")
  rownames(stats.df) <- c("Correlation of scores with factors  ","Multiple R square of scores with factors ","Minimum correlation of factor score estimates")
  print(round(stats.df,digits))
  #cat("\nMeasures of factor score adequacy             ",names(eigenvalues))
  #cat("\nCorrelation of scores with factors            ",round(sqrt(x$stats$R2),digits))
  #cat("\nMultiple R square of scores with factors      " ,round(x$stats$R2,digits))
  #cat("\nMinimum correlation of factor score estimates ", round(2*x$stats$R2 -1,digits),"\n")
   
   cat("\n Total, General and Subset omega for each subset\n")
   colnames(x$omega.group) <- c("Omega total for total scores and subscales","Omega general for total scores and subscales ", "Omega group for total scores and subscales")
   print(round(t(x$omega.group),digits))
   
   #now, see if there are any confidence intervals to report
   if(!is.null(xx$ci)) {
  	 cat("\n Estimates and bootstrapped confidence intervals\n")
  	 li <- data.frame(lower=xx$ci$ci[,1],estimate=xx$ci$means,upper=xx$ci$ci[,2])
  	 li[1,2] <- x$omega_h
  	 li[2,2] <- x$alpha
  	 li[3,2] <- x$omega.tot
  	 li[4,2] <- x$G6
  	 li[5,2] <- x$omega.lim
  	 

   print(li,digits=digits)} 
  }

   
   
   
   
   "print.psych.omegaSem" <- 
   function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 
     if(is.null(cut)) cut <- .2
	 cat( x$title,"\n") 
	 cat("Call: ")
     print(x$Call)
   print.psych.omega(x$omegaSem,digits=digits,all=all,cut=cut,sort=sort,...)
   cat("\n Omega Hierarchical from a confirmatory model using sem = ", round(x$omega.efa$omega,digits))
    cat("\n Omega Total  from a confirmatory model using sem = ", round(x$omega.efa$omega.tot,digits),"\n")
   cat("With loadings of \n")
   loads <- x$omega.efa$cfa.loads
   class(loads) <- NULL
   nfactor <- ncol(loads)
             tn <- colnames(loads)
	        loads <- data.frame(loads)
	        colnames(loads) <- tn  #this seems weird, but otherwise we lose the F* name
	       
           load.2 <- loads     
	        h2 <- round(rowSums(load.2^2),digits)
	        loads <- round(loads,digits)
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")  
         	fx[abs(loads)< cut] <- paste(rep(" ", nc), collapse = "")
            h2 <- round(rowSums(load.2^2),digits)
         	u2 <- 1 - h2
         	
         	p2 <- loads[,ncol(loads)]
         	mp2 <- mean(p2)
         	vp2 <- var(p2)
         	p2 <- round(p2,digits)
	    	print(cbind(fx,h2,u2),quote="FALSE")
	    
	   loads <- as.matrix(load.2) 
	   eigenvalues <- diag(t(loads) %*% loads)
       cat("\nWith eigenvalues of:\n")
       ev.rnd <- round(eigenvalues,digits)
       print(ev.rnd,digits=digits)
   }
