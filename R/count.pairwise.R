#drastically simplified, March 14, 2009 from two loops to 1 matrix operation
#modified July 2, 2013 to allow not counting the diagonal
 "count.pairwise" <-
function (x, y=NULL,diagonal=TRUE) 
{
   if(is.null(y)) {n <- t(!is.na(x)) %*% (!is.na(x)) } else { n <- t(!is.na(x)) %*% (!is.na(y)) } 
   if(!diagonal) diag(n) <- NA
       return(n) }
    
pairwiseDescribe <- function(x,diagonal=FALSE) {
cp <- count.pairwise(x,diagonal=diagonal)
cp <- as.vector(cp[lower.tri(cp)])
describe(cp,skew=FALSE)
}  
    
