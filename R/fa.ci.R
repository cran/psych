#revised 11/5/14 to handle sorted fa data with bootstrapped confidence intervals 
"print.psych.fa.ci" <-
function(x,digits=2,all=FALSE,...)  {
cat("Factor Analysis with confidence intervals using method = ",x$f$fm )
  # cat("\nCall: ")
   print(x$Call)
   class(x) <- c("psych","fa")
   print(x)
  nfactors <-dim(x$loadings)[2]
   c("\n Confidence intervals\n")
   if(is.null(x[["ci"]])) { lc <- lci <- data.frame(unclass(x$loadings),x$cis$ci)} else { lc <- lci <- data.frame(unclass(x$loadings),x$ci)}
  
  
   for(i in 1:nfactors) {
   lci[,(i-1)*3 +2 ] <- lc[,i] 
   lci[,(i-1)*3 +1 ] <- lc[,i+nfactors] 
   lci[,(i-1)*3 +3 ] <- lc[,i+nfactors*2] 
    }
    colnames(lci) <- paste(rep(c("low","coeff","upper"),nfactors),sep="")
   for(i in 1:nfactors) { colnames(lci)[(i-1)*3+2] <- colnames(x$loadings)[i] }
   
   cat("\n Coefficients and bootstrapped confidence intervals \n")
   print (round(lci,digits=digits))
  if(!is.null(x$Phi)) { phis <- x$Phi[lower.tri(x$Phi)]
   cci <- data.frame(lower=x$cis$ci.rot$lower,estimate = phis,upper= x$cis$ci.rot$upper)
    cnR <- abbreviate(colnames(x$loadings),minlength=5) 
      k <- 1
     for(i in 1:(nfactors-1)) {for (j in (i+1):nfactors) {
      rownames(cci)[k] <- paste(cnR[i],cnR[j],sep="-")
      k<- k +1 }}  #added 10/4/14
   cat("\n Interfactor correlations and bootstrapped confidence intervals \n")
   print(cci,digits=digits)}
  
}



