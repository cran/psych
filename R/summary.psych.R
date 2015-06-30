#Modified 1/10/14 to convert to switch 
"summary.psych" <-
function(object,digits=2,items=FALSE,...) { 

#figure what we are trying to summarize
#omega, ICLUST, score.clusters,cluster.cor
#faBy

#if(!is.null(object$title)) { cat("\nSummary of an analysis of ",object$title)}
	
#figure out which psych function called us 

if(length(class(object)) > 1)  { value <- class(object)[2] }

 if(value=="principal") value <- "fa"
 if(value=="score.items") value <- "scores"
 if(value=="cluster.loadings") value <- "cluster.cor"
  if(value=="mat.regress") value <- "mat.reg"
   if(value=="set.cor") value <- "setCor"
    if(value=="mat.reg") value <- "setCor"
 
switch(value, 


 
iclust = { cat("ICLUST (Item Cluster Analysis)") 
 cat("Call: ")
     print(object$call)
    cat( object$title,"\n") 
 	cat("\nPurified Alpha:\n")
	print(object$purified$alpha,digits)
	cat("\n Guttman Lambda6* \n")
	print(object$G6,digits)
	cat("\nOriginal Beta:\n")
	print(object$beta,digits)
	cat("\nCluster size:\n")
	print(object$purified$size,digits)

if(!is.null(object$purified$cor)) {cat("\nPurified scale intercorrelations\n reliabilities on diagonal\n correlations corrected for attenuation above diagonal: \n")
print(object$purified$corrected,digits)  }

} ,

omega =  {
 cat( object$title,"\n") 
 	cat("Alpha:                ",round(object$alpha,digits),"\n") 
 	cat("G.6:                  ",round(object$G6,digits),"\n")
 	cat("Omega Hierarchical:   " ,round(object$omega_h,digits),"\n")
 	cat("Omega H asymptotic:   " ,round(object$omega.lim,digits),"\n")
 	cat("Omega Total           " ,round(object$omega.tot,digits),"\n")
 numfactors <- dim(object$schmid$sl)[2] -3
  eigenvalues <- diag(t(object$schmid$sl[,1:numfactors]) %*% object$schmid$sl[,1:numfactors])
       cat("\nWith eigenvalues of:\n")
       print(eigenvalues,digits=2)
   maxmin <- max(eigenvalues[2:numfactors])/min(eigenvalues[2:numfactors])
   gmax <- eigenvalues[1]/max(eigenvalues[2:numfactors])
  # cat("\ngeneral/max " ,round(gmax,digits),"  max/min =  ",round(maxmin,digits),"\n")
   cat("The degrees of freedom for the model is",object$schmid$dof," and the fit was ",round(object$schmid$objective,digits),"\n")
   	if(!is.na(object$schmid$n.obs)) {cat("The number of observations was ",object$schmid$n.obs, " with Chi Square = ",round(object$schmid$STATISTIC,digits), " with prob < ", round(object$schmid$PVAL,digits),"\n")}
   	
    if(!is.null(object$stats$rms)) {cat("\nThe root mean square of the residuals is ", round(object$stats$rms,digits),"\n") }
     if(!is.null(object$stats$crms)) {cat("The df corrected root mean square of the residuals is ", round(object$stats$crms,digits),"\n") }
    if(!is.null(object$schmid$RMSEA)) {cat("\nRMSEA and the ",object$schmid$RMSEA[4]  ,"confidence intervals are ",round(object$schmid$RMSEA[1:3],digits+1))  }
   	if(!is.null(object$schmid$BIC)) {cat("\nBIC = ",round(object$schmid$BIC,digits))}	
   	  if(!is.null(object$ECV))  cat("Explained Common Variance of the general factor = ", round(object$ECV,digits),"\n")
   	 cat("\n Total, General and Subset omega for each subset\n")
   colnames(object$omega.group) <- c("Omega total for total scores and subscales","Omega general for total scores and subscales ", "Omega group for total scores and subscales")
   print(round(t(object$omega.group),digits)) 	
           }, 

scores =  { #also score.items
cat("Call: ")
print(object$Call)
if(object$raw) {
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (unstandardized) alpha on the diagonal \n corrected correlations above the diagonal:\n") } else {
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (standardized) alpha on the diagonal \n corrected correlations above the diagonal:\n") } 

	 print(object$corrected,digits) 
	 result <- object$corrected
	 },
 
overlap =  { 
cat("Call: ")
print(object$Call)
cat("\nScale intercorrelations adjusted for item overlap")
cat("\nScale intercorrelations corrected for attenuation \n raw correlations (corrected for overlap) below the diagonal, (standardized) alpha on the diagonal \n corrected (for overlap and reliability) correlations above the diagonal:\n") 
	 print(object$corrected,digits) 
	 result <- object$corrected
	 },
	 
vss = {
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
 },
	 
cluster.cor = { 
cat("Call: ")
print(object$Call)
cat("\nScale intercorrelations corrected for attenuation \n raw correlations below the diagonal, (standardized) alpha on the diagonal \n corrected correlations above the diagonal:\n")
	 print(object$corrected,digits) 
     result <- object$corrected

},

fa =  {
  cat("\nFactor analysis with Call: ")
   print(object$Call)
   
 	nfactors <- dim(object$loadings)[2]
    objective <- object$criteria[1]
     if(!is.null(objective)) {    cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.")
    cat("\nThe degrees of freedom for the model is",object$dof," and the objective function was ",round(objective,digits),"\n") 
   	if(!is.na(object$n.obs)) {cat("The number of observations was ",object$n.obs, " with Chi Square = ",round(object$STATISTIC,digits), " with prob < ", signif(object$PVAL,digits),"\n")}
   }


    if(!is.null(object$rms)) {cat("\nThe root mean square of the residuals (RMSA) is ", round(object$rms,digits),"\n") }
    if(!is.null(object$crms)) {cat("The df corrected root mean square of the residuals is ", round(object$crms,digits),"\n") }
    
   
   	if(!is.null(object$TLI)) {cat("\nTucker Lewis Index of factoring reliability = ",round(object$TLI,digits+1))}
  
   	if(!is.null(object$RMSEA)) {cat("\nRMSEA index = ",round(object$RMSEA[1],digits+1), " and the", (1- object$RMSEA[4])*100,"% confidence intervals are ",round(object$RMSEA[2:3],digits+1))  }
   
   	if(!is.null(object$BIC)) {cat("\nBIC = ",round(object$BIC,digits))}
  if(!is.null(object$Phi)) {
   if(object$fn == "principal") {cat ("\n With component correlations of \n" ) } else {cat ("\n With factor correlations of \n" )}
       colnames(object$Phi) <- rownames(object$Phi) <- colnames(object$loadings)
       print(round(object$Phi,digits))}
},

faBy = {cat("\nFactor analysis within groups with Call: ")
 print(object$Call)
 print(object$mean.loading,digits=digits)
 print(object$mean.phi,digits = digits)
 
 
 
 },
      

items= { 
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
		print(object$item.cor,digits) } 

	if (!is.null(object$p.sorted$sorted)) {
 		cat("\nItem by Cluster Structure matrix:\n")
 		print(object$p.sorted$sorted,digits) }
 
 	if (!is.null(object$purified$pattern)) {
 		cat("\nItem by Cluster Pattern matrix:\n")
		 print(object$purified$pattern,digits) }
		 
   if(vss) {
      cat("\nVelicer MAP\n")
      print(object$map,digits)
       cat("\nVery Simple Structure Complexity 1\n")
       print(object$cfit.1,digits)
       cat("\nVery Simple Structure Complexity 2\n")
       print(object$cfit.2,digits)
      }
	
	},  #end if items
	
alpha=  {
cat("\nReliability analysis ",object$title," \n")

print(object$total,digits=digits)
},
 
setCor = {
   if(object$raw) {cat("\nMultiple Regression from raw data \n")} else {
            cat("\nMultiple Regression from matrix input \n")}

              print(object$Call)
            cat("\nMultiple Regression from matrix input \n")
           cat("\nBeta weights \n")
           print(object$beta,digits)
           cat("\nMultiple R \n") 
           print(object$R,digits)
            cat("\nMultiple R2 \n") 
           print(object$R2,digits)
             cat("\nCohen's set correlation R2 \n") 
           print(object$Rset,digits)
           cat("\nSquared Canonical Correlations\n")
           print(object$cancor2,digits)
           
           },
           
irt.fa = {
     cat("\nItem Response Theory using factor analysis with Call: ")
   print(object$Call)
   
 	nfactors <- dim(object$fa$loadings)[2]
    objective <- object$fa$criteria[1]
     if(!is.null(objective)) {    cat("\nTest of the hypothesis that", nfactors, if (nfactors == 1)  "factor is" else "factors are", "sufficient.")
    cat("\nThe degrees of freedom for the model is",object$fa$dof," and the objective function was ",round(objective,digits),"\n") 
   	if(!is.na(object$fa$n.obs)) {cat("The number of observations was ",object$fa$n.obs, " with Chi Square = ",round(object$fa$STATISTIC,digits), " with prob < ", signif(object$fa$PVAL,digits),"\n")}
  
    if(!is.null(object$fa$rms)) {cat("\nThe root mean square of the residuals (RMSA) is ", round(object$fa$rms,digits),"\n") }
    if(!is.null(object$fa$crms)) {cat("The df corrected root mean square of the residuals is ", round(object$fa$crms,digits),"\n") }
    
   
   	if(!is.null(object$fa$TLI)) cat("\nTucker Lewis Index of factoring reliability = ",round(object$fa$TLI,digits+1))}
  
   	if(!is.null(object$fa$RMSEA)) {cat("\nRMSEA index = ",round(object$fa$RMSEA[1],digits+1), " and the", (1- object$fa$RMSEA[4])*100,"% confidence intervals are ",round(object$fa$RMSEA[2:3],digits+1))  }
   
   	if(!is.null(object$fa$BIC)) {cat("\nBIC = ",round(object$fa$BIC,digits))}
    if(!is.null(object$fa$Phi)) {
   if(object$fa$fn == "principal") {cat ("\n With component correlations of \n" ) } else {cat ("\n With factor correlations of \n" )}
       colnames(object$fa$Phi) <- rownames(object$fa$Phi) <- colnames(object$fa$loadings)
       print(round(object$fa$Phi,digits))}
},

describeData = {   cat('n.obs = ', object$n.obs, "of which ", object$complete.cases," are complete cases. Number of variables = ",object$nvar," of which all are numeric is ",object$all.numeric,"\n")
                 
        } 
)  #end of switch
#invisible(result)
   }
  
