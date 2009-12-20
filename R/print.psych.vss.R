"print.psych.vss" <-
function(x,digits=2,all=FALSE,cut=NULL,sort=FALSE,...)  {
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
 cat("\nThe Velicer MAP criterion achieves a minimum of ")
 vss.map <- round(min(x$map) ,digits) 
 cat(vss.map," with " ,which.min(x$map), " factors\n ") 
  cat("\nVelicer MAP\n")
      print(round(x$map,digits))
       cat("\nVery Simple Structure Complexity 1\n")
       print(round(x$cfit.1,digits))
       cat("\nVery Simple Structure Complexity 2\n")
       print(round(x$cfit.2,digits))
   }  #end print.psych.vss
   