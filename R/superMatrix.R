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
 
 
 
 
 #a helper function to take the output of scoreOverlap 
#and form the matrix of (scores by item) x (scores by item)  correlations 
#so that we can then do bestScales on the items related to the score
#Created 12/23/21

"superCor" <- function(x,y=NULL,xy=NULL) {
if(length(class(x)) > 1)  { value <- class(x)[2] } 
if(value %in% c("score.items","overlap")) {
scores <- x
x <- scores$cor
xy <- scores$item.cor}
nx <- NCOL(x)
ny <- NROW(xy)
if(is.null(y)) {nc <- nx
     nr <- nx + ny} else {nc <- nx + NCOL(y)
     nr <- nc
     if(NCOL(y)!= NROW(xy)) stop("Rows of xy matrix do not match rows of y.  Did you remember to specify select=FALSE when using scoreOverlap or scoreItem?  ")}

R <- matrix(NA,nr,nc)
R[1:nx,1:nx] <- x
R[(nx+1):nr,1:nx] <- xy
if(!is.null(y)) {R[1:nx,(nx+1):nc] <- t(xy)
R[(nx+1):nr,(nx+1):nc] <- y}
colnames(R) <- c(colnames(x),colnames(y))
rownames(R ) <- c(colnames(x),rownames(xy))
return(R)}
