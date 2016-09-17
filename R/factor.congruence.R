
#modified October 9, 2015 to add the NA option 
#January 27, 2014  added fa.congruence to clean up calls

#modified March 12 to allow for a list of factor solutions

"factor.congruence" <-
function (x,y=NULL,digits=2,use=NULL,structure=FALSE) {
   fa.congruence(x=x,y=y,digits=digits,use=use,structure=structure) }

"fa.congruence" <-
 function (x,y=NULL,digits=2,use=NULL,structure=FALSE) {
if(is.null(y)&& is.list(x)) {
	n <- length(x)
		for (i in 1:n) {
			xi <- x[[i]]
			if(length(class(xi)) > 1)  { if(class(xi)[2] =='omega')  {xi <- xi$schmid$sl
   												xi<- as.matrix(xi[,1:(ncol(xi)-2)])}}
			if (!(is.matrix(xi) )) {if(!is.null(xi$loadings)) {if(structure) {xi <- xi$Structure} else {xi <- xi$loadings}} else {xi <- as.matrix(xi)}}
if(i==1) {xg <- xi} else {xg <- cbind(xg,xi)} 
}
x <- xg
if(is.null(y)) y <- xg
}  else {
if(length(class(x)) > 1)  { if(class(x)[2] =='omega')  {x <- x$schmid$sl
  		 x <- as.matrix(x[,1:(ncol(x)-2)])}}
if(length(class(y)) > 1)  { if(class(y)[2] =='omega')  {y <- y$schmid$sl
  	 y <- as.matrix(y[,1:(ncol(y)-2)])}}
 
 if (!is.matrix(x)) {if(!is.null(x$loadings)) {if(structure) {x <- x$Structure} else {x <- x$loadings}} else {x <- as.matrix(x)}   }
   if (!is.matrix(y)) {if(!is.null(y$loadings))  {if(structure) {y <- y$Structure} else { y <- y$loadings }} else {y <- as.matrix(y)}}
   }
   
   if(any(is.na(x) | any(is.na(y) ))) {warning("Some loadings were missing.")
        if(!is.null(use)) {message("Analysis is  done on complete cases") 
     if(any(is.na(x))) {
        xc <- x[complete.cases(x),]
        y <- y[complete.cases(x),]
        x <- xc
        }
     if (any(is.na(y))) {
       yc <- y[complete.cases(y),]
       x <- x[complete.cases(y),]
       y <- yc}
     }
     else {warning("Check your data or rerun with the  use = complete option")}
     }
     
   
      
  nx <- dim(x)[2]
  ny <- dim(y)[2]
  cross<- t(y) %*% x   #inner product will have dim of ny * nx
   sumsx<- sqrt(1/diag(t(x)%*%x))   
   sumsy<- sqrt(1/diag(t(y)%*%y)) 

   result<- matrix(rep(0,nx*ny),ncol=nx)
   result<-  round(sumsy * (cross * rep(sumsx, each = ny)),digits)
  
   return(t(result))
   }
