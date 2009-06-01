 "Yule2poly.matrix" <-
function(x,v) {
if (!require(polycor)) {stop("I am sorry, you need to have loaded the polycor package")}
 sizex <- dim(x)[2]
 if (!is.vector(v)) v <- as.vector(v)
 nv <- length(v)
 sizey <- dim(x)[1]
 for (i in 1:sizex) { xc <- v[i]
   for (j in 1:i) {
   	 	yc  <-  v[j]
   	 	if(x[i,j] > .999) {x[i,j] <- 1} else 
  		 {x[i,j] <- Yule2poly(x[i,j],c(xc,yc))}
  		x[j,i] <- x[i,j]
 			 }
 	}
  return(x)  
 }

   
   
   
   
 "phi2poly.matrix" <-
function(x,v) {
if (!require(polycor)) {stop("I am sorry, you need to have loaded the polycor package")}
 sizex <- dim(x)[2]
 if (!is.vector(v)) v <- as.vector(v)
 nv <- length(v)
 sizey <- dim(x)[1]
 for (i in 1:sizex) { xc <- v[i]
   for (j in 1:i) {
   	 	yc  <-  v[j]
   	 	if(x[i,j] > .999) {x[i,j] <- 1} else 
  		 {x[i,j] <- phi2poly(x[i,j],xc,yc)}
  		x[j,i] <- x[i,j]
 			 }
 	}
  return(x)  
 }


"Yule2phi.matrix" <-
function(x,v) {
if (!require(polycor)) {stop("I am sorry, you need to have loaded the polycor package")}
 sizex <- dim(x)[2]
 if (!is.vector(v)) v <- as.vector(v)
 nv <- length(v)
 sizey <- dim(x)[1]
 for (i in 1:sizex) { xc <- v[i]
   for (j in 1:i) {
   	 	yc  <-  v[j]
   	 	if(x[i,j] > .999) {x[i,j] <- 1} else 
  		 {x[i,j] <- Yule2phi(x[i,j],c(xc,yc))}
  		x[j,i] <- x[i,j]
 			 }
 	}
  return(x)  
 }
