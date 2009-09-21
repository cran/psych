"print.psych.fa" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...)  {
if(!is.null(x$fn) ) {if(x$fn == "principal") {cat("Principal Components Analysis") } else {
 cat("Factor Analysis using method = ",x$fm )}}
   cat("\nCall: ")
   print(x$Call)
     
     
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
	    	 fx.1 <- fx[,1,drop=FALSE]    #drop = FALSE  preserves the rownames for single factors
	    	 fx.2 <- fx[,3:(2+ncol)]
	    	 load.2 <- loads[,3:(ncol+2)]
         	fx.2[abs(load.2) < cut] <- paste(rep(" ", nc), collapse = "")
         	fx <- data.frame(V=fx.1,fx.2)
         	if(dim(fx)[2] <3) colnames(fx) <- c("V",colnames(x$loadings)) #for the case of one factor
         	if(nfactors > 1) {if(is.null(x$Phi)) {h2 <- rowSums(load^2)} else {h2 <- diag(load %*% x$Phi %*% t(load)) }} else {h2 <-load.2^2}
         	u2 <- round(1 - h2,digits)
         	h2 <- round(h2,digits)
         	
	    	print(cbind(fx,h2,u2),quote="FALSE")
 		
      	   #adapted from print.loadings
      	  if(is.null(x$Phi)) {if(nfactors > 1)  {vx <- colSums(load.2^2) } else {vx <- sum(load.2^2)
      	                                                     names(vx) <- colnames(x$loadings)
      	                                                     }} else {vx <- diag(x$Phi %*% t(load) %*% load)
      	                                                      names(vx) <- colnames(x$loadings)}
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
    cat("\nThe degrees of freedom for the model is",x$dof," and the objective function was ",round(objective,digits),"\n") 
   	if(!is.na(x$n.obs)) {cat("The number of observations was ",x$n.obs, " with Chi Square = ",round(x$STATISTIC,digits), " with prob < ", signif(x$PVAL,digits),"\n")}

cat("\nFit based upon off diagonal values =", round(x$fit.off,2))
if(!is.null(x$R2)) { stats.df <- t(data.frame(sqrt(x$R2),x$R2,2*x$R2 -1))
 
 rownames(stats.df) <- c("Correlation of scores with factors  ","Multiple R square of scores with factors ","Minimum correlation of factor score estimates")}


 if(!is.null(x$R2)) {cat("\nMeasures of factor score adequacy             \n")
 print(round(stats.df,digits))}

# if(!is.null(x$R2)) {cat("\nMeasures of factor score adequacy             ",colnames(x$loadings)  )
#cat("\nCorrelation of scores with factors           ",round(sqrt(x$R2),digits))
# cat("\nMultiple R square of scores with factors      " ,round(x$R2,digits))
 # cat("\nMinimum correlation of factor score estimates ", round(2*x$R2 -1,digits)) }
 # cat("\nValidity of unit weighted factor scores       ",round(x$valid,digits),"\n")
  
} 
 }  #end of fac.pa
