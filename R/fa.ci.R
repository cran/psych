
 
"print.psych.fa.ci" <-
function(x,digits=2,all=FALSE,...)  {
cat("Factor Analysis with confidence intervals using method = ",x$f$fm )
   cat("\nCall: ")
   print(x$Call)
   print(x$fa)
  nfactors <-dim(x$fa$loadings)[2]
   c("\n Confidence intervals\n")
  
   lc <- lci <- data.frame(unclass(x$fa$loadings),x$ci)
  
   for(i in 1:nfactors) {
   lci[,(i-1)*3 +2 ] <- lc[,i] 
   lci[,(i-1)*3 +1 ] <- lc[,i+nfactors] 
   lci[,(i-1)*3 +3 ] <- lc[,i+nfactors*2] 
    }
    colnames(lci) <- paste(rep(c("low","coeff","upper"),nfactors),sep="")
   for(i in 1:nfactors) { colnames(lci)[(i-1)*3+2] <- colnames(x$fa$loadings)[i] }
   
   cat("\n Coefficients and bootstrapped confidence intervals \n")
   print (round(lci,digits=digits))
  if(!is.null(x$fa$Phi)) { phis <- x$fa$Phi[lower.tri(x$fa$Phi)]
   cci <- data.frame(lower=x$ci.rot$lower,estimate = phis,upper= x$ci.rot$upper)
   cat("\n Interfactor correlations and bootstrapped confidence intervals \n")
   print(cci,digits=digits)}
   
}





 
 
 

 
 
 
