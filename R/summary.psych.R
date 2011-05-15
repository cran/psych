"summary.psych" <-
function(object,digits=2,items=FALSE,...) { 

#figure what we are trying to summarize
#omega, ICLUST, score.clusters,cluster.cor


#if(!is.null(object$title)) { cat("\nSummary of an analysis of ",object$title)}
	
#figure out which psych function called us 
 iclust <- omega <- vss <- scores <- fa <- alpha <- cluster.cor <- mat.reg <- irt.fa <-  FALSE 
 
 if(length(class(object)) > 1)  { 
   if(class(object)[2] =='sim')  sim <- TRUE
   if(class(object)[2] =='vss')  vss <- TRUE
   if(class(object)[2] =='iclust')  iclust <- TRUE
   if(class(object)[2] =='omega')  omega <- TRUE
   if(class(object)[2] =='fa')  fac.pa <- TRUE
   if(class(object)[2] =='principal') fac.pa <- TRUE
   if(class(object)[2] == 'alpha') alpha <- TRUE
   if(class(object)[2] == 'score.items') scores <- TRUE
   if(class(object)[2] == 'cluster.cor') cluster.cor <- TRUE
   if(class(object)[2] == 'cluster.loadings') cluster.cor <- TRUE
   if(class(object)[2] == 'mat.regress') mat.reg <- TRUE
   if(class(object)[2] == 'irt.fa') irt.fa <- TRUE
     } 
result <- NULL
if(!is.null(object$fn)) {fa <- TRUE}
 
 if(vss) {
 if(object$title!="Very Simple Structure") {
 cat("\nVery Simple Structure of ", object$title,"\n") } else {cat("\nVery Simple Structure\n")} 
 cat("VSS complexity 1 achieves a maximimum of ")
 vss.max <- round(max(object$cfit.1) ,digits) 
 cat(vss.max," with " ,which.max(object$cfit.1), " factors\n") 
 cat("VSS complexity 2 achieves a maximimum of ")
  vss.max <- round(max(object$cfit.2) ,digits) 
 cat(vss.max," with " ,which.max(object$cfit.2), " factors\n") 
 cat("\nThe Velicer MAP criterion achieves a minimum of ")
 vss.map <- round(max(object$map) ,digits) 
 cat(vss.map," with " ,which.min(object$map), " factors\n ") 
 }
 
if(iclust) { cat("ICLUST (Item Cluster Analysis)") 
 cat("Call: ")
     print(object$call)
    cat( object$title,"\n") 
 	cat("\nPurified Alpha:\n")
	print(object$purified$alpha)
	cat("\n Guttman Lambda6  * \n")
	print(object$G6,digits)
	cat("\nOriginal Beta:\n")
	print(object$beta)
	cat("\nCluster size:\n")
	print(object$purified$size)

if(!is.null(object$purified$cor)) {cat("\nPurified scale intercorrelations\n reliabilities on diagonal\n correlations corrected for attenuation above diagonal: \n")
print(object$purified$corrected)  }

} 

if(omega) {
 cat( object$title,"\n") 
 cat("Alpha: ",round(object$alpha,digits),"\n") 
 cat("Omega Hierarchical:  " ,round(object$omega_h,digits),"\n")
 cat("Omega Total:  " ,round(object$omega.tot,digits),"\n")
 numfactors <- dim(object$schmid$sl)[2] -2
  eigenvalues <- diag(t(object$schmid$sl[,1:numfactors]) %*% object$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=2)
   maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
   gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
   cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits),"\n")
   cat("The degrees of freedom for the model is",object$schmid$dof," and the fit was ",round(object$schmid$objective,digits),"\n")
   	if(!is.na(object$schmid$n.obs)) {cat("The number of observations was ",object$schmid$n.obs, " with Chi Square = ",round(object$schmid$STATISTIC,digits), " with prob < ", round(object$schmid$PVAL,digits),"\n")}
   	
    if(!is.null(object$stats$rms)) {cat("\nThe root mean square of the residuals is ", round(object$stats$rms,digits),"\n") }
     if(!is.null(object$stats$crms)) {cat("The df corrected root mean square of the residuals is ", round(object$stats$crms,digits),"\n") }
    if(!is.null(object$schmid$RMSEA)) {cat("\nRMSEA and the ",object$schmid$RMSEA[4]  ,"confidence intervals are ",round(object$schmid$RMSEA[1:3],digits+1))  }
   	if(!is.null(object$schmid$BIC)) {cat("\nBIC = ",round(object$schmid$BIC,digits))}	
   	
   	
   	
           } 

if(scores ) { 
cat("Call: ")
print(object$Call)
if(object$raw) {
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (unstandardized) alpha on the diagonal \n corrected correlations above the diagonal:\n") } else {
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (standardized) alpha on the diagonal \n corrected correlations above the diagonal:\n") } 

	 print(object$corrected,digits) 
	 result <- object$corrected
	 }
 
if(cluster.cor ) { 
cat("Call: ")
print(object$Call)
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (standardized) alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(object$corrected,digits) 
     result <- object$corrected

}

if(fa) {
  cat("\nFactor analysis with Call: ")
   print(object$Call)
   
 	nfactors <- dim(object$loadings)[2]
    objective <- object$criteria[1]
     if(!is.null(objective)) {    cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.")
    cat("\nThe degrees of freedom for the model is",object$dof," and the objective function was ",round(objective,digits),"\n") 
   	if(!is.na(object$n.obs)) {cat("The number of observations was ",object$n.obs, " with Chi Square = ",round(object$STATISTIC,digits), " with prob < ", signif(object$PVAL,digits),"\n")}
   }


    if(!is.null(object$rms)) {cat("\nThe root mean square of the residuals is ", round(object$rms,digits),"\n") }
    if(!is.null(object$crms)) {cat("The df corrected root mean square of the residuals is ", round(object$crms,digits),"\n") }
    
   
   	if(!is.null(object$TLI)) cat("\nTucker Lewis Index of factoring reliability = ",round(object$TLI,digits+1))}
  
   	if(!is.null(object$RMSEA)) {cat("\nRMSEA index = ",round(object$RMSEA[1],digits+1), " and the", (1- object$RMSEA[4])*100,"% confidence intervals are ",round(object$RMSEA[2:3],digits+1))  }
   
   	if(!is.null(object$BIC)) {cat("\nBIC = ",round(object$BIC,digits))

}


if(items) { 
    if(omega) {
       cat("\nSchmid Leiman Factor loadings:\n")
       print(object$schmid$sl)
       numfactors <- dim(object$schmid$sl)[2] -2
       eigenvalues <- diag(t(object$schmid$sl[,1:numfactors]) %*% object$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=digits)
       
       }
    
	if(!is.null(object$item.cor) ) {
		cat("\nItem by scale correlations:\n")
		print(object$item.cor) } 

	if (!is.null(object$p.sorted$sorted)) {
 		cat("\nItem by Cluster Structure matrix:\n")
 		print(object$p.sorted$sorted) }
 
 	if (!is.null(object$purified$pattern)) {
 		cat("\nItem by Cluster Pattern matrix:\n")
		 print(object$purified$pattern) }
		 
   if(vss) {
      cat("\nVelicer MAP\n")
      print(round(object$map,2))
       cat("\nVery Simple Structure Complexity 1\n")
       print(round(object$cfit.1,2))
       cat("\nVery Simple Structure Complexity 2\n")
       print(round(object$cfit.2,2))
      }
	
	}  #end if items
	
if(alpha) {
cat("\nReliability analysis ",object$title," \n")

print(object$total,digits=digits)
}
 
  if(mat.reg) { cat("\nMultiple Regression for matrix input \nCall: ")
              print(object$Call)
            cat("\nMultiple Regression from matrix input \n")
           cat("\nBeta weights \n")
           print(round(object$beta,digits))
           cat("\nMultiple R \n") 
           print(round(object$R,digits))
            cat("\nMultiple R2 \n") 
           print(round(object$R2,digits))
             cat("\nCohen's set correlation R2 \n") 
           print(round(object$Rset,digits))
           
           }
           
  if(irt.fa) {
   cat("Item Response Analysis using Factor Analysis = ")
   cat("\nCall: ")
   print(object$Call)
  print(round(object$coefficients,digits))
  print(object$stats,digits)
  }
invisible(result)
   }
  
