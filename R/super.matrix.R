 "super.matrix" <- 
 function(x,y) {
if(is.vector(x)) {x <- matrix(x)
                  colnames(x) <- "X"
                  if(dim(x)[1] <2) {rownames(x) <- "X"}
                  } else {if (is.null(colnames(x))) colnames(x) <- paste("X",1:dim(x)[2],sep="")}
if(is.vector(y)) {y <- matrix(y)
                  colnames(y) <- "Y"
                  if(dim(y)[1]<2) {rownames(y) <- "Y"}
                  } else {if (is.null(colnames(y))) colnames(y) <- paste("Y",1:dim(x)[2],sep="")}
 fillx <- rbind(x,matrix(0,ncol=dim(x)[2],nrow=dim(y)[1]))
 filly <- rbind(matrix(0,ncol=dim(y)[2],nrow=dim(x)[1]),y) 
 xy <- cbind(fillx,filly)
 colnames(xy) <- c(colnames(x),colnames(y))
 rownames(xy) <- c(rownames(x),rownames(y))
 return(xy)}
 