 "superMatrix" <- 
 function(x,y=NULL) {if(is.list(x)) {
                 if(is.null(y)) { y <- x[-1]} else { y <- list(x[-1],y) }
                  x <- x[[1]]
                  }  
                if(is.list(y)) {
	            if (length(y) > 1)  {
		         x <- superMatrix(x,y[[1]])
                 xy <- superMatrix(x,y[-1])
              } else {y <- y[[1]]
                 xy <- superMatrix(x,y)}
                 } else {
if(is.vector(x)) {x <- matrix(x)
                  colnames(x) <- "X"
                  if(dim(x)[1] <2) {rownames(x) <- "X"} else {rownames(x) <- paste("Vx",1:dim(x)[1],sep="") }
                  } else {if (is.null(colnames(x))) colnames(x) <- paste("X",1:dim(x)[2],sep="")
                          if (is.null(rownames(x))) rownames(x) <- paste("Vx",1:dim(x)[1],sep="")}
if(is.vector(y)) {y <- matrix(y)
                  colnames(y) <- "Y"
                  if(dim(y)[1]<2) {rownames(y) <- "Y"} else {rownames(y) <- paste("Vy",1:dim(y)[1],sep="") }
                  } else {if (is.null(colnames(y))) colnames(y) <- paste("Y",1:dim(y)[2],sep="")
                          if (is.null(rownames(y))) rownames(y) <- paste("Vy",1:dim(y)[1],sep="")}
 fillx <- rbind(x,matrix(0,ncol=dim(x)[2],nrow=dim(y)[1]))
 filly <- rbind(matrix(0,ncol=dim(y)[2],nrow=dim(x)[1]),y) 
 xy <- cbind(fillx,filly)
 colnames(xy) <- c(colnames(x),colnames(y))
 rownames(xy) <- c(rownames(x),rownames(y))
 }
 
 return(xy)}
 #fixed June 21, 2009 to add rownames of matrices if necessary
 #modified June 8, 2012 to add list input option
 
  "super.matrix" <- 
 function(x,y) {
 .Deprecated("super.matrix", msg = "super.matrix is deprecated.  Please use the superMatrix function")
xy <- superMatrix(x,y)
 return(xy)}
 #fixed June 21, 2009 to add rownames of matrices if necessary
 #modified June 8, 2012 to add list input option