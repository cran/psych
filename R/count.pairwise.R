#drastically simplified, March 14, 2009 from two loops to 1 matrix operation
 "count.pairwise" <-
function (x, y=NULL) 
{
   if(is.null(y)) {n <- t(!is.na(x)) %*% (!is.na(x)) } else { n <- t(!is.na(x)) %*% (!is.na(y)) } 
    return(n) }
    
    
