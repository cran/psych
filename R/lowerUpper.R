#written June 6, 2012
#note that lower.tri and upper.tri return the matrix in a different order
"lowerUpper" <-
function(lower,upper=NULL,diff=FALSE) {
if(is.null(upper)) {upper <- lower  #return two from one
    upper[lower.tri(upper)] <- t(upper)[lower.tri(t(upper))]
    lower <- t(lower)
    lower[lower.tri(lower)] <- t(lower)[lower.tri(lower)]
    result <- list(lower=lower,upper=upper)
    } else {
    if(nrow(lower) !=ncol(lower)) {stop("lower matrix must be square")}
if(nrow(upper) !=ncol(upper)) {stop("upper matrix must be square")}
if(nrow(lower) !=ncol(upper)) {stop("lower and upper matrices must have the same dimensions")}
result <- lower
colnames(result) <- colnames(upper)
rownames(result) <-rownames(lower)
if(diff) upper <- lower - upper
result [lower.tri(result)] <- upper[lower.tri(upper)] 
result <- t(result)
diag(result) <- NA}
return(result)}
#revised Oct 6, 2013 to pick up row names and column names from the two matrices
