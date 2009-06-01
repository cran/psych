"print.psych.iclust" <- 
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

   cat("ICLUST (Item Cluster Analysis)")
   cat("\nCall: ")
     print(x$call)
     
     
   
 	cat("\nPurified Alpha:\n")
	print(x$purified$alpha,digits)
		cat("\nG6* reliability:\n")
	print(x$G6,digits)
	cat("\nOriginal Beta:\n")
	print(x$beta)
	cat("\nCluster size:\n")
	print(x$purified$size,digits)

 
if(sort) { cat("\nItem by Cluster Structure matrix: Sorted by loading \n")
 load <- x$sorted$sorted
  if(is.null(cut)) cut <- .3
	    ncol <- dim(load)[2]-3
	    fx <- as.matrix(format(load,digits=digits))
	     nc <- nchar(fx[1,4], type = "c")
	     fx.1 <- fx[,1:3]
	     fx.2 <- fx[,4:(3+ncol)]
	     load.2 <- load[,4:(ncol+3)]
         fx.2[abs(load.2)< cut] <- paste(rep(" ", nc), collapse = "")
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
         fx[abs(load) < cut] <- paste(rep(" ", nc), collapse = "")
	    print(fx,quote="FALSE")
 		#print(unclass(x$loadings)) 
 		eigenvalues <- diag(t(x$loadings) %*% x$loadings)
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
 		}
 if(!is.null(x$purified$cor)) {cat("\nPurified scale intercorrelations\n reliabilities on diagonal\n correlations corrected for attenuation above diagonal: \n")
print(x$purified$corrected)  }

   }# end of if ICLUST