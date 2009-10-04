"summary.psych" <-
function(object,digits=2,items=FALSE,...) { 

#figure what we are trying to summarize
#omega, ICLUST, score.clusters,cluster.cor


#if(!is.null(object$title)) { cat("\nSummary of an analysis of ",object$title)}
	
#figure out which psych function called us 
 iclust <- omega <- vss <- scores <- fa <- alpha <- cluster.cor <- FALSE 
 
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
     } 

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
           } 

if(scores ) { 
cat("Call: ")
print(object$Call)
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (unstandardized) alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(object$corrected,digits) 
	 }
 
if(cluster.cor ) { 
cat("Call: ")
print(object$Call)
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (standardized) alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(object$corrected,digits) 
 

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
	
if(alpha) {
cat("\nReliability analysis ",object$title," \n")

print(object$total,digits=digits)
}
 
}