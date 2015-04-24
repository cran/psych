"print.psych.fa" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,suppress.warnings=TRUE,...)  {
if(!is.matrix(x) && !is.null(x$fa) && is.list(x$fa)) x <-x$fa   #handles the output from fa.poly
if(!is.null(x$fn) ) {if(x$fn == "principal") {cat("Principal Components Analysis") } else {
 cat("Factor Analysis using method = ",x$fm )}}
   cat("\nCall: ")
   print(x$Call)
     
   	load <- x$loadings  
 
 if(is.null(cut)) cut <- 0   #caving into recommendations to print all loadings
 
    #but, if we are print factors of covariance matrices, they might be very small
     #  cut <- min(cut,max(abs(load))/2)   #removed following a request by  Reinhold Hatzinger
     
 	nitems <- dim(load)[1]
 	nfactors <- dim(load)[2]
 		if(sum(x$uniqueness) + sum(x$communality) >  nitems) {covar <- TRUE} else {covar <- FALSE}
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
                if(max(abs(load) > 1.0) && !covar) cat('\n Warning: A Heywood case was detected. \n')
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
            colnames(fx) <- colnames(x$loadings)}
         	if(nfactors > 1) {if(is.null(x$Phi)) {h2 <- rowSums(load.2^2)} else {h2 <- diag(load.2 %*% x$Phi %*% t(load.2)) }} else {h2 <-load.2^2}
         	if(!is.null(x$uniquenesses)) {u2 <- x$uniquenesses[u2.order]}  else {u2 <- (1 - h2)}
         	#h2 <- round(h2,digits)
         	vtotal <- sum(h2 + u2)
           if(isTRUE(all.equal(vtotal,nitems))) {
           cat("Standardized loadings (pattern matrix) based upon correlation matrix\n")
           com <- x$complexity[u2.order] # u2.order added 9/4/14
           if(!is.null(com)) { print(cbind(fx,h2,u2,com),quote="FALSE",digits=digits)} else {
            print(cbind(fx,h2,u2),quote="FALSE",digits=digits) } } else {
	    	cat("Unstandardized loadings (pattern matrix) based upon covariance matrix\n") 
	    	print(cbind(fx,h2,u2,H2=h2/(h2+u2),U2=u2/(h2+u2)),quote="FALSE",digits=digits)}
	   
 		
      	   #adapted from print.loadings
      	  if(is.null(x$Phi)) {if(nfactors > 1)  {vx <- colSums(load.2^2) } else {vx <- sum(load.2^2)
      	                                             }} else {vx <- diag(x$Phi %*% t(load) %*% load)
      	                                                      }
          
          names(vx) <- colnames(x$loadings)
          varex <- rbind("SS loadings" =   vx)
          varex <- rbind(varex, "Proportion Var" =  vx/vtotal)
         if (nfactors > 1) {
                      varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/vtotal))
                       varex <- rbind(varex, "Proportion Explained"=  vx/sum(vx))
                      varex <- rbind(varex, "Cumulative Proportion"=  cumsum(vx/sum(vx))) 
                      }
                      
             cat("\n")
            print(round(varex, digits))
           
                      
          #now, if we did covariances show the standardized coefficients as well
	    	
	    	if(!isTRUE(all.equal(vtotal,nitems))) {  #total variance accounted for is not just the number of items in the matrix
	    	cat('\n Standardized loadings (pattern matrix)\n')
	    	 
	    	 fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")
	    	 fx.1 <- fx[,1,drop=FALSE]    #drop = FALSE  preserves the rownames for single factors
	    	 fx.2 <- round(loads[,3:(2+ncol)]/sqrt(h2+u2),digits)
	    	 load.2 <- loads[,3:(ncol+2)]/sqrt(h2+u2)
	    	 
	    	
         	fx.2[abs(load.2) < cut] <- paste(rep(" ", nc), collapse = "")
         	fx <- data.frame(V=fx.1,fx.2)
         	if(dim(fx)[2] <3) colnames(fx) <- c("V",colnames(x$loadings)) #for the case of one factor
         	if(nfactors > 1) { h2 <-h2/(h2+u2)} else {h2 <-h2/(h2+u2)}
         	u2 <- (1 - h2)
         	
	    	print(cbind(fx,h2,u2),quote="FALSE",digits=digits)


  if(is.null(x$Phi)) {if(nfactors > 1)  {vx <- colSums(load.2^2) } else {vx <- diag(t(load) %*% load)
                                                                         vx <- vx*nitems/vtotal 
      	                                                         }} else {vx <- diag(x$Phi %*% t(load) %*% load)
      	                                                           vx <- vx*nitems/vtotal }
      	  names(vx) <- colnames(x$loadings)
          
           
          varex <- rbind("SS loadings" =   vx)
          varex <- rbind(varex, "Proportion Var" =  vx/nitems)
           if (nfactors > 1) {varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/nitems))
                              varex <- rbind(varex, "Cum. factor Var"=  cumsum(vx/sum(vx)))}
    cat("\n") 
    print(round(varex, digits))
     
    	    	}
  
    
    if(!is.null(x$Phi))  { 
      if(!is.null(x$fn) ) { if(x$fn == "principal") {cat ("\n With component correlations of \n" ) } else {cat ("\n With factor correlations of \n" )}}
       colnames(x$Phi) <- rownames(x$Phi) <- colnames(x$loadings)
       print(round(x$Phi,digits))} else {
       if(!is.null(x$rotmat)) {
             U <- x$rotmat
             ui <- solve(U)
           Phi <- t(ui) %*% ui
          Phi <- cov2cor(Phi) 
         if(!is.null(x$fn) ) { if(x$fn == "principal") {cat ("\n With component correlations of \n" ) } else {cat ("\n With factor correlations of \n" )}}
       colnames(Phi) <- rownames(Phi) <- colnames(x$loadings)
       print(round(Phi,digits))
            } }
       
       if(!is.null(x$complexity)) cat("\nMean item complexity = ",round(mean(x$complexity),1))     
       objective <- x$criteria[1]
    if(!is.null(objective)) { if(!is.null(x$fn) ) { if(x$fn == "principal") {  cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "component is" else "components are", "sufficient.\n")} else { cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.\n")}}
  if(x$fn != "principal") {
    if(!is.null(x$null.dof)) {cat("\nThe degrees of freedom for the null model are ",x$null.dof, " and the objective function was ",round(x$null.model,digits),...)}
    if(!is.null(x$null.chisq)) {cat(" with Chi Square of " ,round(x$null.chisq,digits)) }
    cat("\nThe degrees of freedom for the model are",x$dof," and the objective function was ",round(objective,digits),"\n",...) 
     }
     
    if(!is.null(x$rms)) {cat("\nThe root mean square of the residuals (RMSR) is ", round(x$rms,digits),"\n") }
    if(!is.null(x$crms)) {cat("The df corrected root mean square of the residuals is ", round(x$crms,digits),"\n",...) }
    
     if((!is.null(x$nh)) && (!is.na(x$nh))) {cat("\nThe harmonic number of observations is " ,round(x$nh)) }
     if((!is.null(x$chi)) && (!is.na(x$chi))) {cat(" with the empirical chi square ", round(x$chi,digits), " with prob < ", signif(x$EPVAL,digits),"\n" ,...)  }
   	 if(x$fn != "principal") { 
   	 if(!is.na(x$n.obs)) {cat("The total number of observations was ",x$n.obs, " with MLE Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits),"\n",...)}
  
     
   	if(!is.null(x$TLI)) cat("\nTucker Lewis Index of factoring reliability = ",round(x$TLI,digits+1))}
   	if(!is.null(x$RMSEA)) {cat("\nRMSEA index = ",round(x$RMSEA[1],digits+1), " and the", (1- x$RMSEA[4])*100,"% confidence intervals are ",round(x$RMSEA[2:3],digits+1),...)  }
   	if(!is.null(x$BIC)) {cat("\nBIC = ",round(x$BIC,digits))}
}
 	if(!is.null(x$fit)) cat("\nFit based upon off diagonal values =", round(x$fit.off,digits))

if ((!is.null(x$fn)) && (x$fn != "principal")) {
if(!is.null(x$R2)) { stats.df <- t(data.frame(sqrt(x$R2),x$R2,2*x$R2 -1))

 rownames(stats.df) <- c("Correlation of scores with factors  ","Multiple R square of scores with factors ","Minimum correlation of possible factor scores ")
         colnames(stats.df) <- colnames(x$loadings)
 } else {stats.df <- NULL}
 badFlag <- FALSE
 #however, if the solution is degenerate, don't print them
  

   if( (is.null(x$R2)) || (any(max(x$R2,na.rm=TRUE) > (1 + .Machine$double.eps)) )) {badFlag <- TRUE
      if (!suppress.warnings) {
	cat("\n WARNING, the factor score fit indices suggest that the solution is degenerate. Try a different method of factor extraction.\n")
	 warning("the factor score fit indices suggest that the solution is degenerate\n")}
	 } else {
	
	 if(!is.null(stats.df)) { cat("\nMeasures of factor score adequacy             \n")
      print(round(stats.df,digits))}

	if(!is.null(x$stats$R2) && !badFlag) {cat("\nMeasures of factor score adequacy             ",colnames(x$loadings)  )
	cat("\nCorrelation of scores with factors           ",round(sqrt(x$R2),digits))
	cat("\nMultiple R square of scores with factors      " ,round(x$R2,digits))
	  cat("\nMinimum correlation of factor score estimates ", round(2*x$R2 -1,digits)) }
	  }
 }
 result <- list(Vaccounted=varex)
 invisible(result) 
} 
   #end of print.psych.fa

#modified November 22, 2010 to get the communalities correct for sorted loadings, but does this work for covariances?
#modified November 18, 2012 to print the empirical chi squares
#modified October 13, 2013 to add the invisibile return of varex.