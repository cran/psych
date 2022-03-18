#modified June 25, 2018 to handle omegaSem output as well
#modified October 9, 2015 to add the NA option 
#January 27, 2014  added fa.congruence to clean up calls

#modified March 12 to allow for a list of factor solutions
#Modified December 11, 2019 to use inherits rather than class 
"factor.congruence" <-
function (x,y=NULL,digits=2,use=NULL,structure=FALSE) {
   fa.congruence(x=x,y=y,digits=digits,use=use,structure=structure) }

"fa.congruence" <-
 function (x,y=NULL,digits=2,use=NULL,structure=FALSE) {
 direct <- extend <-  esem <- factanal <- other  <- NA
  obnames <- cs(fa, omega, omegaSem, directSl, direct, omegaDirect, principal, iclust,extend,esem, factanal)

if(is.null(y) && is.list(x)) {
	n <- length(x)
		for (i in 1:n) {
			xi <- x[[i]]
			if(length(class(xi)) > 1)  {
			   cln <- inherits(xi, obnames, which=TRUE)
			   if (any(cln > 1)) {cln <- obnames[which(cln >0)]} else {cln <- "other"}} else {cln <- "other"}

			    switch(cln,
			      fa = {if(structure) {xi <- xi$Structure} else {xi <- xi$loadings}},
			      omega = {xi <- xi$schmid$sl
   				 		  xi <- as.matrix(xi[,1:(ncol(xi)-2)])},
   				 omegaSem = {xi <- xi$omega.efa$cfa.loads},
   				 directSl = {xi <- xi$direct},
   				 direct = {xi <- xi$direct},
   				  omegaDirect = {xi <- xi$loadings},
   				 principal = {xi <- xi$loadings},
   				 iclust = {xi <- xi$loadings},
   				  extend = {xi <- xi$loadings},
   				  esem = {xi <- xi$loadings},
   				 other = {if(inherits(xi, "factanal")) {xi <- xi$loadings} else {xi <- as.matrix(xi)}}
   				 )
   				 if(i==1) {xg <- xi} else {xg <- cbind(xg,xi)} 
   				}		  
		x <- xg
		
if(is.null(y)) y <- xg
}  else {
if(length(class(x)) > 1) {#cln <- class(x)[2]} else {cln <- "other"}
               cln <- inherits(x, obnames, which=TRUE)
			   if (any(cln > 1)) {cln <- obnames[which(cln >0)]} else {cln <- "other"}} else {cln <- "other"}  #fixe March 3, 2020
			   
			    switch(cln,
			     fa = {if(structure) {x <- x$Structure} else {x <- x$loadings}},
			      omega = {x <- x$schmid$sl
   				 		  x <- as.matrix(x[,1:(ncol(x)-2)])},
   				 omegaSem = {x <- x$omega.efa$cfa.loads},
   				 directSl = {x <- x$direct},
   				 direct = {x <- x$direct},
   				 omegaDirect = {x <- x$loadings},
   				 principal = {x <- x$loadings},
   				 iclust = {x <- x$loadings},
   				  extend = {x <- x$loadings},
   				  esem = {x <- x$loadings},
   				  other = {if(inherits(x,  "factanal")) {x <- x$loadings} else {x <- as.matrix(x)}}
   				 )
   				}		  
  		 
if(length(class(y)) > 1) {   #{ cln <- class(y)[2] } else {cln <- "other"}
              cln <- inherits(y, obnames, which=TRUE)
			   if (any(cln > 1)) {cln <- obnames[which(cln >0)]} else {cln <- "other"}
			   } else {cln <- "other"}
			    switch(cln,
			       fa = {if(structure) {y <- y$Structure} else {y <- y$loadings}},
			      omega = {y <- y$schmid$sl
   				 		  y <- as.matrix(y[,1:(ncol(y)-2)])},
   				 omegaSem = {y <- y$omega.efa$cfa.loads},
   				 directSl = {y <- y$direct},
   				 direct = {y <- y$direct},
   				  omegaDirect = {y <- y$loadings},
   				  principal = {y <- y$loadings},
   				  esem = {y <- y$loadings},
   				    extend = {y <- y$loadings},
   				 iclust = {y <- y$loadings},
   				   other = {if(inherits(y, "factanal")) {y <- y$loadings} else {y <- as.matrix(y)}}
   				 )

   
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
   
   
   
   
 #find the generalized congruence coefficient 
 #normalized cross products
 #if zero centered data, this is the correlation
 #if centered on the response midpoint, this is the Cohen cs 
 
#handles missing data
#find the congruence coefficient
#handles missing data
congruence <- function(x,y=NULL){
if(is.null(y)) y <- x
nvarx <- NCOL(x)
nvary <- NCOL(y)
if(nvarx == 1) x <- as.matrix(x)
if(nvary ==1) y <-as.matrix(y)

C <- matrix(nrow=nvary,ncol=nvarx)
Cx <- rep(0,nvarx)
Cy <- rep(0,nvary)
for(i in 1:nvarx) {
 for(j in 1:nvary) {
 C[j,i]  <- sum(x[,i,drop=FALSE] * y[,j,drop=FALSE],na.rm=TRUE)
  Cx[i] <- sum(x[,i,drop=FALSE]^2,na.rm=TRUE)}}
 
 for(j in 1:nvary){
   Cy[j] <- sum(y[,j,drop=FALSE]^2,na.rm=TRUE)}
 
     
 if(nvarx> 1) {dcx <- sqrt(diag(1/Cx))} else {dcx <- sqrt(1/Cx)}
 if(nvary > 1) {dcy <- sqrt(diag(1/Cy))} else {dcy <- sqrt(1/Cy)}

 C <- dcy %*% C %*% dcx
 rownames(C) <- colnames(y)
 colnames(C) <- colnames(x)
 class(C) <- c("psych","congruence")
 return(C)
}

cohen.profile <- function(x,y=NULL,M=NULL) {
if(is.null(y)) y <- x
if(is.null(M)) {
min.scale<- min(x,y,na.rm=TRUE)
max.scale <- max(x,y,na.rm=TRUE)
M <- (max.scale + min.scale )/2}
congruence(x=x-M,y=y-M)}

distance <- function(x,y=NULL,r=2) {
if(is.null(y)) y <- x
nvarx <- NCOL(x)
nvary <- NCOL(y)
C <- matrix(nrow=nvary,ncol=nvarx)
Cx <- rep(0,nvarx)
Cy <- rep(0,nvary)
for(i in 1:nvarx) {
 for(j in 1:nvary) {
 if(r==1) {C[j,i] <- sum(abs(x[,i,drop=FALSE] - y[,j,drop=FALSE])) } else {
 C[j,i]  <- exp(log(sum((x[,i,drop=FALSE] - y[,j,drop=FALSE])^r,na.rm=TRUE))/r)}
  }
  }
 rownames(C) <- colnames(y)
 colnames(C) <- colnames(x)
  C
  }
