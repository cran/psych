 "super.matrix" <- 
 function(x,y) {
 fillx <- rbind(x,matrix(0,ncol=dim(x)[2],nrow=dim(y)[1]))
 filly <- rbind(matrix(0,ncol=dim(y)[2],nrow=dim(x)[1]),y) 
 xy <- cbind(fillx,filly)
 colnames(xy) <- c(colnames(x),colnames(y))
 rownames(xy) <- c(rownames(x),rownames(y))
 return(xy)}
 