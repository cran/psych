"print.psych.stats" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...) { 

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
#end of print.psych.stats


   
   
 