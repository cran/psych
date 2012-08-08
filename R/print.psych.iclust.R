"print.psych.iclust" <- 
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

   cat("ICLUST (Item Cluster Analysis)")
   cat("\nCall: ")
     print(x$call)
     
    if((!is.null(x$purify)) && x$purify) {
   
 	cat("\nPurified Alpha:\n")
	print(x$purified$alpha,digits)
		cat("\nG6* reliability:\n")
	print(x$purified$G6,digits)
	cat("\nOriginal Beta:\n")
	print(x$beta,digits)
	cat("\nCluster size:\n")
	print(x$purified$size,digits)
} else 
  {
  	cat("\noriginal Alpha:\n")
	print(x$alpha,digits)
		cat("\nG6* reliability:\n")
	print(x$G6,digits)
	cat("\nOriginal Beta:\n")
	print(x$beta,digits)
	cat("\nCluster size:\n")
	print(x$size,digits)}
if(sort) { cat("\nItem by Cluster Structure matrix: Sorted by loading \n")
 load <- x$sorted$sorted
  if(is.null(cut)) cut <- .3
	    ncol <- dim(load)[2]-3
	    load[4:(ncol+3)] <- round(load[4:(ncol+3)],digits)
	    fx <- as.matrix(format(load,digits=digits))
	     nc <- nchar(fx[1,4], type = "c")
	     fx.1 <- fx[,1:3]
	     fx.2 <- format(fx[,4:(3+ncol)],digits)
	     load.2 <- load[,4:(ncol+3)]
         fx.2[abs(load.2)< cut] <- paste(rep(" ", nc), collapse = "")
         fx <- data.frame(fx.1,fx.2)
	    print(fx,quote="FALSE")
 		
 		eigenvalues <- diag(t(x$pattern) %*% x$loadings)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
 		}
 	else {
	 cat("\nItem by Cluster Structure matrix:\n")
	    load <- unclass(x$loadings)
	    load <- round(load,digits)
	    fx <- format(load,digits=digits)
	     nc <- nchar(fx[1,1], type = "c")
         fx[abs(load) < cut] <- paste(rep(" ", nc), collapse = "")
        if(is.matrix(x$clusters)) {
         clust <- colnames(x$clusters)[apply(abs(x$clusters),1,which.max)]
         pclust <- colnames(x$p.sorted$clusters)[apply(abs(x$p.sorted$clusters),1,which.max)]
           clust.fx <- data.frame(O=clust,P=pclust,fx)} else {clust.fx <- fx}
         
	    print(clust.fx,quote="FALSE")
 		#print(unclass(x$loadings)) 
 		eigenvalues <- diag(t(x$pattern) %*% x$loadings)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
 		}
 if(!is.null(x$purified$cor)) {cat("\nPurified scale intercorrelations\n reliabilities on diagonal\n correlations corrected for attenuation above diagonal: \n")
print(round(x$purified$corrected,digits=digits))  }

cat("\nCluster fit = ",round(x$fit$clusterfit,digits), "  Pattern fit = ",round(x$fit$patternfit,digits), " RMSR = ",round(x$fit$patternrmse,digits), "\n")

   }# end of if ICLUST