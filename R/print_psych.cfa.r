#Developed 12/21/2025
"print_psych.cfa" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,suppress.warnings=TRUE,...)  {
# if(!is.null(x$irt)) x<-x$irt
if(!is.null(x$ci)) {ci <- x$ci
 x <- x$cf
 }
 
if(!is.matrix(x) && !is.null(x$fa) && is.list(x$fa)) x <-x$fa   #handles the output from fa.poly
if(!is.matrix(x))  {if (!is.null(x$fn) ) {if(x$fn == "principal") {cat("Principal Components Analysis") } else {
 cat("Factor Analysis using method = ",x$fm )}  }} else {load <- x
                                                         class(x)<- NULL
                                                         print(round(x,digits))
                                                         
                                                         return()} 
   cat("\nCall: ")
   print(x$Call)
    
   	rload <- load <- x$loadings  
   	 if(is.null(cut)) cut <- 0   #caving into recommendations to print all loadings
 
   	nitems <- dim(load)[1]
 	nfactors <- dim(load)[2]
   	ncol <- dim(load)[2]
   	   cat('Confirmatory factor (Structure)  loadings\n')
   
     rloads <- round(cbind(load, x$communality),digits)
     colnames(rloads)[nfactors+1] <- "h2"
	    	fx <- format(rloads,digits=digits)     	
	
	print(rloads)
	if(nfactors >1 ) {
	     fx <- format(rloads,digits=digits)
	    	#nc <- nchar(rloads[1,3], type = "c")    #dropped because we don't have communalities
	    	 fx.1 <- fx[,1,drop=FALSE]    #drop = FALSE  preserves the rownames for single factors
	    	 #fx.2 <- fx[,3:(2+ncol),drop=FALSE]
	    	 load.2 <- as.matrix(load)
	    	 
	    	 
	    	         	#fx[abs(load[,]) < cut,] <- paste(rep(" ", nc), collapse = "")
	} else {load.2<- as.matrix(load)}
	
	h2 <- x$communality
	#if(nfactors > 1) {if(is.null(x$Phi)) {h2 <- rowSums(load.2^2)} else {h2 <- diag(load.2 %*% x$Phi %*% t(load.2)) }} else {h2 <-load.2^2}
         #	if(!is.null(x$uniquenesses)) {u2 <- x$uniquenesses[u2.order]}  else 
         	{u2 <- (1 - h2)}
vtotal <- sum(h2 + u2)	

if(nfactors > 1)  {vx <- colSums(load^2) } else {vx <- sum(load^2)}
 names(vx) <- colnames(load)
          varex <- rbind("SS loadings" =   vx)
          varex <- rbind(varex, "Proportion Var" =  vx/vtotal)
         if (nfactors > 1) {
                      varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/vtotal))
                       varex <- rbind(varex, "Proportion Explained"=  vx/sum(vx))
                      varex <- rbind(varex, "Cumulative Proportion"=  cumsum(vx/sum(vx))) 
                      }     	                                            	    	
	    	print(round(varex, digits))
	    	
	    if(!is.null(x$Phi))	{
	    	 cat('With correlations of \n')
	    	lowerMat(x$Phi)
	    	}
	    	
	    	#now print the results fron stats   (if they exist,  they do not for omegaStats)
	    	if(!is.null(x$stats)) {
	    	xx <- x$stats  #allowing calls using the calls from print.psych.fa
	if(!is.null(x$dof))  {cat("\nThe degrees of freedom of the model are ",x$dof)}
 if(!is.null(xx$rms)) {cat("\nThe root mean square of the residuals (RMSR) is ", round(xx$rms,digits),"\n") }
    if(!is.null(xx$crms)) {cat("The df corrected root mean square of the residuals is ", round(xx$crms,digits),"\n",...) }
    
     if((!is.null(xx$nh)) && (!is.na(xx$nh))) {cat("\nThe harmonic n.obs is " ,round(xx$nh)) }
     if((!is.null(xx$chi)) && (!is.na(xx$chi))) {cat(" with the empirical chi square ", round(xx$chi,digits), " with prob < ", signif(xx$EPVAL,digits),"\n" ,...)  }
   	  
   	 if(!is.na(xx$n.obs)) {cat("The total n.obs was ",xx$n.obs, " with Likelihood Chi Square = ",round(xx$STATISTIC,digits), " with prob < ", signif(xx$PVAL,digits),"\n",...)}
  
     
   	if(!is.null(xx$TLI)) {cat("\nTucker Lewis Index of factoring reliability = ",round(xx$TLI,digits+1))}
   	if(!is.null(xx$RMSEA)) {cat("\nRMSEA index = ",round(xx$RMSEA[1],digits+1), " and the", (xx$RMSEA[4])*100,"% confidence intervals are ",round(xx$RMSEA[2:3],digits+1),...)  }
   	if(!is.null(xx$BIC)) {cat("\nBIC = ",round(xx$BIC,digits))}	    	
}
#get information from factor stats

if(!is.null(x$stats$R2)) { stats.df <- t(data.frame(sqrt(x$stats$R2),x$stats$R2,2*x$stats$R2 -1)) 

 rownames(stats.df) <- c("Correlation of (regression) scores with factors  ","Multiple R square of scores with factors ","Minimum correlation of possible factor scores ")
       #  rownames(stats.df) <- colnames(xx$loadings)
} else  {stats.df <- NULL}

if(!is.null(x$omegaStats)) {    #the normal case
if(is.null(x$X)) {X <- ci$X} else {X<-x$X}
nfactors <- ncol(X)
fnames <- colnames(X)

 if(!is.null(stats.df)) { cat("\nMeasures of factor score adequacy          \n")
      print(round(stats.df,digits))}
   

if(is.null(x$gr.ci)) {  #this is the case of CFA output
cat("\n Loadings and their 95% confidence intervals\n")
for(i in 1:nfactors){
cat("\n Factor ",i, fnames[i], "\n")
print(round(ci$ci[X[,i]>0,],digits))
}
} else { #the case of bifactor output
cat("\n general factor loadings and their 95% confidence intervals\n")
print(round(x$gf.ci,digits))
cat("\n group factor loadings and their 95% confidence intervals\n")

for(i in 1:nfactors){
cat("\n Factor ",i, fnames[i], "\n")
print(round(x$gr.ci[X[,i]>0,],digits))
}


}  
}   

#omega stats
 #two cases:  results from  CFA.bifactor or from omegaStats
     
  if(!is.null(x$omegaStats)) {cat("\nMeasures of internal consistency based on the bifactor model:  \n")
           cat("    omega_bi    = ", round(x$omegaStats$omega_h,digits)," Note that this is not omega_hierarchical nor an estimate of the general factor.", "\n    alpha       = " ,
           round(x$omegaStats$alpha,digits), "\n    omega_total = ",round(x$omegaStats$omega_total,digits))
          cat("\n Total, General and Subset omega for each factor\n") 
           colnames(x$omegaStats$omega.group) <- c("Omega total for total scores and subscales","Omega general for total scores and subscales ", "Omega group for total scores and subscales")
   print(round(t(x$omegaStats$omega.group),digits))
           }  else {if(!is.null(x$omega_h)) {
           cat("\nMeasures of internal consistency based on the bifactor model:  \n")
           cat("    omega_bi    = ", round(x$omega_h,digits)," Note that this is not omega_hierarchical nor an estimate of the general factor.", "\n    alpha       = " ,
           round(x$alpha,digits), "\n    omega_total = ",round(x$omega_total,digits))
          cat("\n Total, General and Subset omega for each factor\n") 
           colnames(x$omega.group) <- c("Omega total for total scores and subscales","Omega general for total scores and subscales ", "Omega group for total scores and subscales")
   print(round(t(x$omega.group),digits))
           }
           
           }
}
