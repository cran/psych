"factor.rotate" <-
function(f,angle,col1=1,col2=2,plot=FALSE,...)  {
   #hand rotate two factors from a loading matrix
   #see the GPArotation package for much more elegant procedures
     if (!is.matrix(f) ) {f <-f$loadings} 
     nvar<- dim(f)[2]
     if(!is.matrix(f)) {if(!is.data.frame(f)) {stop("f must be either a data frame or a matrix")} else {f <- as.matrix(f)} }
     rot<- diag(1,nvar,nvar)
     
     theta<- pi*angle/180
     rot[col1,col1]<-  cos(theta)
     rot[col2,col2]<-  cos(theta)
     rot[col1,col2]<- -sin(theta)
     rot[col2,col1]<-  sin(theta)
     result <- f %*% rot
     if(plot) {fa.plot(result,...)
        abline(a=0,b=tan(-theta),lty="dashed")
        abline(a=0,b=tan(-theta+ pi/2),lty="dashed") }
     return(result) }

