#reorganized, January 18, 2009 to make some what clearer
#should eventually just break into separate print functions
"print.psych" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

iclust <- omega <- vss <- scores <- fac.pa <- principal <- gutt <- sim <- alpha <- describe <- corr.test <- r.test <- cortest <-  cluster.cor <- cluster.loadings <- thurstone <- stats <- FALSE
#first, figure out which psych function was called
if(length(class(x)) > 1)  {
   if(class(x)[2] =='sim')  sim <- TRUE
   if(class(x)[2] =='vss')  vss <- TRUE
   if(class(x)[2] =='iclust')  iclust <- TRUE
   if(class(x)[2] =='omega')  omega <- TRUE
   if(class(x)[2] =='fa')  fac.pa <- TRUE
   if(class(x)[2] =='principal') principal <- fac.pa <- TRUE
   if(class(x)[2] == 'alpha') alpha <- TRUE
   if(class(x)[2] == 'describe') describe <- TRUE
   if(class(x)[2] == "corr.test") corr.test <- TRUE
   if(class(x)[2] == "r.test") r.test <- TRUE
   if(class(x)[2] == "cortest") cortest <- TRUE
   if(class(x)[2] == "score.items") scores <- TRUE
   if(class(x)[2] == "cluster.cor") cluster.cor <- TRUE
   if(class(x)[2] == "cluster.loadings") cluster.loadings <- TRUE
   if(class(x)[2] == "guttman") gutt <- TRUE
   if(class(x)[2] == "thurstone") thurstone <- TRUE
   if(class(x)[2] == "stats") stats <- TRUE
     } 
else {      #old and clunky way of figuring out which function called,   replace as above
   
#these next test for non-psych functions that may be printed using print.psych
if(!is.null(x$communality.iterations)) {fac.pa <- TRUE}
if(!is.null(x$uniquenesses)) {fac.pa <- TRUE}
if(!is.null(x$rotmat)) {fac.pa <- TRUE}
if(!is.null(x$Th)) {fac.pa <- TRUE}

}


if(describe) {if  (length(dim(x))==1) {class(x) <- "list"
              attr(x,"call") <- NULL
              print(x,digits)
                  } else  {class(x) <- "data.frame" 
            print(round(x,digits)) }
      }
         
          
if(corr.test) {cat("Call:")
              print(x$Call)
              cat("Correlation matrix \n")
              print(round(x$r,digits))
              cat("Sample Size \n")
              print(x$n)
              cat("Probability value \n")
             print(round(x$p,digits))
     }

if(r.test) {cat("Correlation tests \n")
            cat("Call:")
              print(x$Call)
              cat( x$Test,"\n")
              if(!is.null(x$t)) {cat(" t value" ,round(x$t,digits),"   with probability <", signif(x$p,digits) )}
              if(!is.null(x$z)) {cat(" z value" ,round(x$z,digits),"   with probability ",  round(x$p,digits) )}               
              if(!is.null(x$ci)) {cat("\n and confidence interval ",round(x$ci,digits) ) }
         }
         
 if(cortest) {cat("Tests of correlation matrices \n")
            cat("Call:")
            print(x$Call)
            cat(" Chi Square value" ,round(x$chi,digits)," with df = ",x$df, "  with probability <", signif(x$p,digits) )
         }
                                
  
if(all) {class(x) <- "list"
         print(x) } else {    #find out which function created the data and then print accordingly

## 
 if(omega) {
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
	   if(sort) {
	      ord <- sort(abs(loads[,1]),decreasing=TRUE,index.return=TRUE)
	      loads[,] <- loads[ord$ix,]
	      loads <- cbind(v=ord$ix,loads)
	   } #end sort 
	    	fx <- format(loads,digits=digits)
	    	nc <- nchar(fx[1,3], type = "c")  
         	fx[abs(loads)< cut] <- paste(rep(" ", nc), collapse = "")
	    	print(fx,quote="FALSE")
	   
       numfactors <- dim(x$schmid$sl)[2] -2
       eigenvalues <- diag(t(x$schmid$sl[,1:numfactors]) %*% x$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
      
  	 maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
  	 gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
   	cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits))
   	cat("\nThe degrees of freedom for the model is",x$schmid$dof," and the fit was ",round(x$schmid$objective,digits),"\n")
   	if(!is.na(x$schmid$n.obs)) {cat("The number of observations was ",x$schmid$n.obs, " with Chi Square = ",round(x$schmid$STATISTIC,digits), " with prob < ", signif(x$schmid$PVAL,digits),"\n")}
   	}

#### VSS outpt
if(vss) {
 if(x$title!="Very Simple Structure") {
 
 cat("\nVery Simple Structure of ", x$title,"\n") } else {cat("\nVery Simple Structure\n")} 
 cat("Call: ")
 print(x$call)
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
    cat("Call: ")
    print(x$Call)
	cat("\n(Unstandardized) Alpha:\n")
	print(x$alpha,digits=digits)
	
  	cat("\nAverage item correlation:\n")
  	print(x$av.r,digits=digits)
	
	cat("\n Guttman 6* reliability: \n")
	print(x$G6,digits=digits)     
	
	 if(iclust) {cat("\nOriginal Beta:\n")
	 print(x$beta,digits) }	 
	           
 	 # cat("\nScale intercorrelations:\n")
	 # print(x$cor,digits=digits)
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	 
	  if(!is.null(x$item.cor) ) {
	   cat("\nItem by scale correlations:\n corrected for item overlap and scale reliability\n" )
	 print(round(x$item.corrected,digits=digits)) } 
  }
  
if(cluster.cor) {
    cat("Call: ")
    print(x$Call)
	cat("\n(Standardized) Alpha:\n")
	print(x$alpha,digits)

  	cat("\nAverage item correlation:\n")
	print(x$av.r,digits)
	cat("\nNumber of items:\n")
	print(x$size)
	
	      
 	 # cat("\nScale intercorrelations:\n")
	 # print(x$cor,digits=digits)
	 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	  }
	
	
if(cluster.loadings) {
    cat("Call: ")
    print(x$Call)
	cat("\n(Standardized) Alpha:\n")
	print(x$alpha,digits)
		cat("\n(Standardized) G6*:\n")
	print(x$G6,digits)
  	cat("\nAverage item correlation:\n")
	print(x$av.r,digits)
	cat("\nNumber of items:\n")
	print(x$size)
 cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(x$corrected,digits) 
	 
   cat("\nItem by scale intercorrelations\n corrected for item overlap and scale reliability\n")
	 print(x$loadings,digits) 
  #cat("\nItem by scale Pattern matrix\n")
	# print(x$pattern,digits) 
	
  }
####
if(alpha) {
cat("\nReliability analysis ",x$title," \n")
cat("Call: ")
print(x$call)
cat("\n ")
print(x$total,digits=digits)
cat("\n Reliability if an item is dropped:\n")
    print(x$alpha.drop,digits=digits)
cat("\n Item statistics \n")
   print(x$item.stats,digits=digits)
}




#### factor.pa output
if(fac.pa) {
 if(is.null(cut)) cut <- .3
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
  #now sort the loadings that have their highest loading on each cluster
  
  		items <- table(loads$cluster)   #how many items are in each cluster?
  		first <- 1
  		item <- loads$item
		for (i in 1:length(items)) {
		if(items[i] > 0 ) {
				last <- first + items[i]- 1
				ord <- sort(abs(loads[first:last,i+2]),decreasing=TRUE,index.return=TRUE)
   				loads[first:last,3:(nfactors+2)] <- load[item[ord$ix+first-1],]
   				loads[first:last,1] <- item[ord$ix+first-1]
   				rownames(loads)[first:last] <- rownames(loads)[ord$ix+first-1]
   		 		first <- first + items[i]  }
          		 }  
         }    #end of sort 		 
    #they are now sorted, don't print the small loadings 
          	ncol <- dim(loads)[2]-2
	    	fx <- format(round(loads,digits=digits))
	    	nc <- nchar(fx[1,3], type = "c")
	    	 fx.1 <- fx[,1]
	    	 fx.2 <- fx[,3:(2+ncol)]
	    	 load.2 <- loads[,3:(ncol+2)]
         	fx.2[abs(load.2)< cut] <- paste(rep(" ", nc), collapse = "")
         	fx <- data.frame(V=fx.1,fx.2)
         	if(dim(fx)[2] <3) colnames(fx) <- c("V",colnames(x$loadings)) #for the case of one factor
	    	print(fx,quote="FALSE")
 		
      	   #adapted from print.loadings
      	  if(nfactors > 1)  {vx <- colSums(load.2^2) } else {vx <- sum(load.2^2)
      	                                                     names(vx) <- colnames(x$loadings)
      	                                                     }
           varex <- rbind("SS loadings" =   vx)
           varex <- rbind(varex, "Proportion Var" =  vx/nitems)
           if (nfactors > 1) 
                      varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/nitems))
  
    cat("\n")
  
    print(round(varex, digits))
    
    if(!is.null(x$Phi))  { 
       cat ("\n With factor correlations of \n" )
       colnames(x$Phi) <- rownames(x$Phi) <- colnames(x$loadings)
       print(round(x$Phi,digits))} else {
       if(!is.null(x$rotmat)) {
             U <- x$rotmat
             ui <- solve(U)
           Phi <- t(ui) %*% ui
          Phi <- cov2cor(Phi) 
           cat ("\n With factor correlations of \n" )
       colnames(Phi) <- rownames(Phi) <- colnames(x$loadings)
       print(round(Phi,digits))
            } }
            
       objective <- x$criteria[1]
     if(!is.null(objective)) {    cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.\n")
    cat("\nThe degrees of freedom for the model is",x$dof," and the fit was ",round(objective,digits),"\n") 
   	if(!is.na(x$n.obs)) {cat("The number of observations was ",x$n.obs, " with Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits),"\n")}
  if(!principal) { cat("\nMeasures of factor score adequacy             ",colnames(x$loadings)  )
  cat("\n Correlation of scores with factors           ",round(sqrt(x$R2),digits))
  cat("\nMultiple R square of scores with factors      " ,round(x$R2,digits))
  cat("\nMinimum correlation of factor score estimates ", round(2*x$R2 -1,digits)) }
  cat("\nValidity of unit weighted factor scores       ",round(x$valid,digits),"\n")

}
          
	 }  #end of fac.pa
 
 ####  ICLUST output
if(iclust) {
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
}
   }# end of if ICLUST


###

 if(gutt) {
  cat("Call: ")
    print(x$Call)
 cat ("\nAlternative estimates of reliability\n")
 cat("Beta = ", x$beta, " This is an estimate of the worst split half reliability")  
 cat ("\nGuttman bounds \nL1 = ",x$lambda.1, "\nL2 = ", x$lambda.2, "\nL3 (alpha) = ", x$lambda.3,"\nL4 (max) = " ,x$lambda.4, "\nL5 = ", x$lambda.5, "\nL6 (smc) = " ,x$lambda.6, "\n")
 cat("TenBerge bounds \nmu0 = ",x$tenberge$mu.0, "mu1 = ", x$tenberge$mu1, "mu2 = " ,x$tenberge$mu2, "mu3 = ",x$tenberge$mu3 , "\n")
 cat("\nalpha of first PC = ", x$alpha.pc, "\nestimated greatest lower bound = ", x$lambda.4,"\n")
 cat("\nbeta estimated by first and second PC = ", round(x$beta.pc,digits), " This is an exploratory statistic \n")
 } 
  
  
 ## 
 if(thurstone) {
cat("Thurstonian scale (case 5) scale values ")
cat("\nCall: ")
     print(x$Call)
print(x$scale)
cat("\n Goodness of fit of model  ", round(x$GF,digits))
 }
  
 if(sim) { if(is.matrix(x)) {x <-unclass(x) 
             round(x,digits) } else {
              cat("Call: ")
              print(x$call)
             cat("\n $model (Population correlation matrix) \n")
             print(x$model,digits)
             if(!is.null(x$reliability)) { cat("\n$reliability (population reliability) \n")
                print(x$reliability,digits) } 
             if(!is.null(x$N) && !is.null(x$r)) {
             cat("\n$r  (Sample correlation matrix  for sample size = ",x$N,")\n")
             print(x$r,digits)}} 
             
      }  #end of the not list condition
  
 if(stats) {
  cat("Call: ")
  print(x$Call)
  nfactors <- x$factors
   objective <- x$criteria[1]
     if(!is.null(objective)) {    cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.\n")
    cat("\nThe degrees of freedom for the model is",x$dof," and the fit was ",round(objective,digits),"\n") 
   	if(!is.na(x$n.obs)) {cat("The number of observations was ",x$n.obs, " with Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits),"\n")}
}
  cat("\nMeasures of factor score adequacy             ",colnames(x$loadings)  )
  cat("\n Correlation of scores with factors           ",round(sqrt(x$R2),digits))
  cat("\nMultiple R square of scores with factors      " ,round(x$R2,digits))
  cat("\nMinimum correlation of factor score estimates ", round(2*x$R2 -1,digits),"\n")
 }

}  #end function
