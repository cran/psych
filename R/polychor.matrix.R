

 "Yule2poly.matrix" <-
function(x,v) {   
.Deprecated("Yule2tetra", msg="This function has been replaced by Yule2tetra")
}

 "phi2poly.matrix" <-
function(x,v) {  
.Deprecated("phi2tetra",msg="This function has been replaced by phi2tetra")
}

"Yule2phi.matrix" <-
function(x,v) {
.Deprecated("Yule2phi",msg="This function has been replaced by Yule2phi")
}
   
   


#revised August 29, 2010 
"poly.mat" <- 
function(x,short=TRUE,std.err=FALSE,ML=FALSE) {
.Deprecated("polychoric",msg="poly.mat is deprecated.  Please use the polychoric function instead.")
return(polychoric(x))	
}

