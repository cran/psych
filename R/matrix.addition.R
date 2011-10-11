#matrix.addition 
"%+%" <- function(x,y) { 
 if(!is.matrix(x)) {
	if(is.vector(x)) {x <- as.matrix(x)} else stop("x must be either a vector or a matrix")}
if(!is.matrix(y)) {
	if(is.vector(y)) {y <- as.matrix(y)} else stop("y must be either a vector or a matrix")}
n.x <- dim(x)[1]
n.y <- dim(y)[2]
n.k <- dim(x)[2]
if (n.k != dim(y)[1]) {warning("Matrices should be comparable")}
#first find sum vectors
x <- rowSums(x,na.rm=FALSE)
y <- colSums(y,na.rm=FALSE)
one <- as.vector(rep(1,n.y))  #to duplicate x n.y times
one.y <- as.vector(rep(1,n.x)) #to duplicate y n.x times
xy  <- x %*% t(one) + t(y %*% t(one.y) )  #sum the vectors in a rectangular array
  return(xy) }
  
 "tr" <- 
function(m) {
if (!is.matrix(m) |(dim(m)[1] != dim(m)[2]) ) stop ("m must be a square matrix")
return(sum(diag(m))) }