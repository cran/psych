"summary.psych" <-
function(object,digits=2,items=FALSE,...) { 

#figure what we are trying to summarize
#omega, ICLUST, score.clusters,cluster.cor


#if(!is.null(object$title)) { cat("\nSummary of an analysis of ",object$title)}
	
#figure out which psych function called us 
 iclust <- omega <- vss <- scores <- fa <-  FALSE 
if (!is.null(object$map)) { vss <- TRUE} 
if(!is.null(object$purified$alpha)) {  iclust <- TRUE} 
if(!is.null(object$omega_h)) {omega <- TRUE}
if(!is.null(object$av.r)) {scores <- TRUE}
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
 
if(iclust) { 
    cat( object$title,"\n") 
 	cat("\nPurified Alpha:\n")
	print(object$purified$alpha)
	cat("\nOriginal Beta:\n")
	print(object$beta)
	cat("\nCluster size:\n")
	print(object$purified$size)

if(!is.null(object$purified$cor)) {cat("\nPurified scale intercorrelations:\n")
print(object$purified$cor)  }

} 

if(omega) {
 cat( object$title,"\n") 
 cat("Alpha: ",object$alpha,"\n") 
 cat("Omega Hierarchical:  " ,object$omega_h,"\n")
 cat("Omega Total:  " ,object$omega.tot,"\n")
 numfactors <- dim(object$schmid$sl)[2] -2
  eigenvalues <- diag(t(object$schmid$sl[,1:numfactors]) %*% object$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=2)
   maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
   gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
   cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits),"\n")
   cat("The degrees of freedom for the model is",object$schmid$dof," and the fit was ",round(object$schmid$objective,digits),"\n")
   	if(!is.na(object$schmid$n.obs)) {cat("The number of observations was ",object$schmid$n.obs, " with Chi Square = ",round(object$schmid$STATISTIC,digits), " with prob < ", round(object$schmid$PVAL,digits),"\n")}
           } 

if(scores & !iclust) { cat("\nAlpha:\n")
print(object$alpha)
 cat("\nScale intercorrelations:\n")
 print(object$cor)  
 

}

if(fa) {print(object$loadings)}


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
 
}