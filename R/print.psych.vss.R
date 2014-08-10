"print.psych.vss" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...)  {
if(!is.null(x$cfit.1)) {
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
 cat("\nThe Velicer MAP achieves a minimum of ")
 vss.map <- round(min(x$map) ,digits) 
 cat(vss.map," with " ,which.min(x$map), " factors ") 
 bic.min <- round(min(x$vss.stats[,"BIC"]),digits)
 cat("\nBIC achieves a minimum of ", bic.min," with ", which.min(x$vss.stats[,"BIC"])," factors\n")
 sabic.min <- round(min(x$vss.stats[,"SABIC"]),digits)
 cat("Sample Size adjusted BIC achieves a minimum of ", sabic.min," with ", which.min(x$vss.stats[,"SABIC"])," factors\n")
 # cat("\nVelicer MAP\n")
  #    print(round(x$map,digits))
  #     cat("\nVery Simple Structure Complexity 1\n")
  #     print(round(x$cfit.1,digits))
  #     cat("\nVery Simple Structure Complexity 2\n")
  #     print(round(x$cfit.2,digits))
 
 temp <- data.frame(vss1=x$cfit.1,vss2=x$cfit.2,map=x$map,x$vss.stats[,1:13])
 cat("\nStatistics by number of factors \n")
 print(temp,digits=digits) } else { if(x$title!="Number of Factors") {
 
 cat("\nNumber of factors of ", x$title,"\n") } else {cat("\nNumber of factors\n")} 
 cat("Call: ")
 print(x$call)
 cat("VSS complexity 1 achieves a maximimum of ")
 vss.max <- round(max(x$vss.stats$cfit.1) ,digits) 
 cat(vss.max," with " ,which.max(x$vss.stats$cfit.1), " factors\n") 
 cat("VSS complexity 2 achieves a maximimum of ")
  vss.max <- round(max(x$vss.stats$cfit.2) ,digits) 
 cat(vss.max," with " ,which.max(x$vss.stats$cfit.2), " factors\n") 
 cat("The Velicer MAP achieves a minimum of ")
 vss.map <- round(min(x$map,na.rm=TRUE) ,digits) 
 cat(vss.map," with " ,which.min(x$map), " factors ") 
 bic.min <- round(min(x$vss.stats[["eBIC"]],na.rm=TRUE),digits)
 cat("\nEmpirical BIC achieves a minimum of ", bic.min," with ", which.min(x$vss.stats[["eBIC"]])," factors\n")
 sabic.min <- round(min(x$vss.stats[["SABIC"]],na.rm=TRUE),digits)
 cat("Sample Size adjusted BIC achieves a minimum of ", sabic.min," with ", which.min(x$vss.stats[["SABIC"]])," factors\n")
 # cat("\nVelicer MAP\n")
  #    print(round(x$map,digits))
  #     cat("\nVery Simple Structure Complexity 1\n")
  #     print(round(x$cfit.1,digits))
  #     cat("\nVery Simple Structure Complexity 2\n")
  #     print(round(x$cfit.2,digits))
 
 temp <- data.frame(vss1=x$vss.stats$cfit.1,vss2=x$vss.stats$cfit.2,map=x$map,x$vss.stats[,1:13])
 cat("\nStatistics by number of factors \n")
 print(temp,digits=digits)}
   }  #end print.psych.vss
   
