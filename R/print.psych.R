"print.psych" <-
function(x,digits=2,all=FALSE,cutoff=0,sort=FALSE,...) { 

iclust <- omega <- vss <- scores <- fac.pa <- gutt <- FALSE

if (!is.null(x$map)) { vss <- TRUE}     
if(!is.null(x$purified$alpha)) { iclust <- TRUE} 
if(!is.null(x$omega_h)) {omega <- TRUE}
if(!is.null(x$av.r)) {scores <- TRUE}
if(!is.null(x$communality.iterations)) {fac.pa <- TRUE}
if(!is.null(x$lambda.1) ) {gutt <- TRUE}

 
if(!is.null(x$scores)) {
 x <- x[-c(1:2)]}
  
if(all) {class(x) <- "list"
         print(x) } else {
 
 if(omega) {
	 cat( x$title,"\n") 
 	cat("Alpha: ",x$alpha,"\n") 
 	cat("Lambda.6:  ",x$lambda.6,"\n")
 	cat("Omega Hierarchical:  " ,x$omega_h,"\n")
 	cat("Omega Total  " ,x$omega.tot,"\n")
            
	cat("\nSchmid Leiman Factor loadings:\n")
       print(x$schmid$sl)
       numfactors <- dim(x$schmid$sl)[2] -2
       eigenvalues <- diag(t(x$schmid$sl[,1:numfactors]) %*% x$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
      
  	 maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
  	 gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
   	cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits),"\n")
   }

if(vss) {
 if(x$title!="Very Simple Structure") {
 cat("\nVery Simple Structure of ", x$title,"\n") } else {cat("\nVery Simple Structure\n")} 
 cat("VSS complexity 1 achieves a maximimum of ")
 vss.max <- round(max(x$cfit.1) ,digits) 
 cat(vss.max," with " ,which.max(x$cfit.1), " factors\n") 
 cat("VSS complexity 2 achieves a maximimum of ")
  vss.max <- round(max(x$cfit.2) ,digits) 
 cat(vss.max," with " ,which.max(x$cfit.2), " factors\n") 
 cat("\nThe Velicer MAP criterion achieves a minimum of ")
 vss.map <- round(max(x$map) ,digits) 
 cat(vss.map," with " ,which.min(x$map), " factors\n ") 
  cat("\nVelicer MAP\n")
      print(round(x$map,digits))
       cat("\nVery Simple Structure Complexity 1\n")
       print(round(x$cfit.1,digits))
       cat("\nVery Simple Structure Complexity 2\n")
       print(round(x$cfit.2,digits))
 }

if(scores) {
	 cat("\nAlpha:\n")
	print(x$alpha)
  	cat("\nAverage item correlation:\n")
	print(x$av.r,digits)
	
if(iclust) {cat("\nOriginal Beta:\n")
		print(x$beta,digits) }	
 	cat("\nScale intercorrelations:\n")
	print(x$cor,digits=digits)
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	 }
	 

if(!is.null(x$item.cor) ) {
	cat("\nItem by scale correlations:\n")
		print(x$item.cor) } 

#if(!is.null(x$sorted)) { load <- x$sorted
if(FALSE) {
	    ncol <-dim(load)[2]-4
	    fx <- format(load,digits=digits)
	     nc <- nchar(fx[1,4], type = "c")
	     fx.1 <- fx[,1:3]
	     fx.2 <- fx[,4:(4+ncol)]
	     load.2 <- load[,4:(ncol+4)]
         fx.2[abs(load.2)< cutoff] <- paste(rep(" ", nc), collapse = "")
         fx <- data.frame(fx.1,fx.2)
	    print(fx,quote="FALSE")
 	   loadings <- load[,4:(4+ncol)]
 	   loadings <- as.matrix(loadings)
 		eigenvalues <- diag(t(loadings) %*% loadings)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)}



if(fac.pa) {
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
          	ncol <-dim(loads)[2]-2
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")
	    	 fx.1 <- fx[,1]
	    	 fx.2 <- fx[,3:(2+ncol)]
	    	 load.2 <- loads[,3:(ncol+2)]
         	fx.2[abs(load.2)< cutoff] <- paste(rep(" ", nc), collapse = "")
         	fx <- data.frame(item=fx.1,fx.2)
	    	print(fx,quote="FALSE")
 		
 			eigenvalues <- diag(t(x$loadings) %*% x$loadings)
       	cat("\nWith eigenvalues of:\n")
      	 print(eigenvalues,digits=digits)  
           
	 } else {   #we don't want to sort
	 print(x,cutoff=cutoff) } 
	 }  #end of fac.pa
 
if(iclust) {
if(sort) { cat("\nItem by Cluster Structure matrix: Sorted by loading \n")
	    load <- x$sorted$sorted
	    ncol <-dim(load)[2]-4
	    fx <- format(load,digits=digits)
	     nc <- nchar(fx[1,4], type = "c")
	     fx.1 <- fx[,1:3]
	     fx.2 <- fx[,4:(4+ncol)]
	     load.2 <- load[,4:(ncol+4)]
         fx.2[abs(load.2)< cutoff] <- paste(rep(" ", nc), collapse = "")
         fx <- data.frame(fx.1,fx.2)
	    print(fx,quote="FALSE")
 		
 		eigenvalues <- diag(t(x$loadings) %*% x$loadings)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
 		}
 	else {
	 cat("\nItem by Cluster Structure matrix:\n")
	    load <- unclass(x$loadings)
	    fx <- format(load,digits=digits)
	     nc <- nchar(fx[1,1], type = "c")
         fx[abs(load) < cutoff] <- paste(rep(" ", nc), collapse = "")
	    print(fx,quote="FALSE")
 		#print(unclass(x$loadings)) 
 		eigenvalues <- diag(t(x$loadings) %*% x$loadings)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
 		}
 
 if (!is.null(x$pattern))  {
    if (sort) {cat("\nItem by Cluster Pattern matrix: sort by loadings is not yet available\n Redo without sort=TRUE to get pattern matrix\n") } else {
       
        cat("\nItem by Cluster Pattern matrix:\n")
		 load <- unclass(x$pattern)
		  fx <- format(load,digits=digits)
	     nc <- nchar(fx[1,1], type = "c")
         fx[abs(load) < cutoff] <- paste(rep(" ", nc), collapse = "")
	    print(fx,quote="FALSE")
		# print(x$pattern) 
		 	eigenvalues <- diag(t(x$pattern) %*% x$pattern)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)}
       }
  
   } #end if !is.null(x$loadings)
   }# end of if ICLUST
   
 if(gutt) {
 cat ("Alternative estimates of reliability \n")
 cat ("Guttman bounds \n L1 = ",x$lambda.1, "  L2 = ", x$lambda.2, "L3 (alpha) = ", x$lambda.3, "L4 (max) = " ,x$lambda.4, "L5 = ", x$lambda.5, "L6 (smc) = " ,x$lambda.6, "\n")
 cat("TenBerge bounds \n mu0 = ",x$tenberge$mu.0, "mu1 = ", x$tenberge$mu1, "mu2 = " ,x$tenberge$mu2, "mu3 = ",x$tenberge$mu3 , "\n")
 cat("Beta = ", x$beta, "        alpha of first PC = ", x$alpha.pc, "estimated glb = ", x$lambda.4,"\n")
 } 
   

}  #end function
